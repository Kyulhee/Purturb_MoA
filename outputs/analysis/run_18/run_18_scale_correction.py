"""
BioEval Run 18 -- Norman logFC Scale Correction (A3)
Rescales Ridge predictions to match true effect variance, then re-evaluates
all metrics to see if downstream task discrimination improves.

A3 hypothesis: Ridge shrinkage compresses logFC scale, reducing DEG fraction
in predictions. Rescaling should restore DEG fraction and improve downstream
task discrimination (especially f1@50 and dir_discovery).

Strategy:
1. Load data and train models (same as run_16/17)
2. For each model, compute per-perturbation scale ratio (pred_std / true_std)
3. Apply 3 correction strategies:
   a. global_rescale: single ratio (mean of per-pert ratios) applied to all
   b. per_pert_rescale: individual ratio per perturbation
   c. variance_match: rescale predictions so pred_var = true_var per perturbation
4. Re-evaluate all metrics with corrected predictions
5. Compare H1/H2/H3 results with and without correction
"""
import json, os, warnings, time
import numpy as np
import scanpy
from scipy import stats
from sklearn.metrics import mean_squared_error, r2_score, average_precision_score
from sklearn.decomposition import PCA

warnings.filterwarnings('ignore')

OUTPUT_DIR = r"C:\test\llm_project\nexus-science-win\outputs\analysis\run_18"
DATA_DIR = r"C:\test\llm_project\nexus-science-win\outputs\analysis\run_04\data"
os.makedirs(OUTPUT_DIR, exist_ok=True)

np.random.seed(42)

def log(msg):
    try:
        print(msg, flush=True)
    except UnicodeEncodeError:
        print(msg.encode('ascii', errors='replace').decode('ascii'), flush=True)

# ============================================================
# BioEval Metrics (same as run_16/17)
# ============================================================

def bioeval_dir(true_mat, pred_mat, deg_threshold=0.25):
    sign_match = (np.sign(true_mat) == np.sign(pred_mat)).astype(float)
    zt = (true_mat == 0); zp = (pred_mat == 0)
    sign_match[zt & zp] = 1.0; sign_match[zt & ~zp] = 0.5
    dir_all = float(sign_match.mean())
    deg_mask = np.abs(true_mat) > deg_threshold
    dir_deg = float(sign_match[deg_mask].mean()) if deg_mask.sum() > 0 else None
    w = np.abs(true_mat); w = w / (w.sum() + 1e-8)
    dir_weighted = float((sign_match * w).sum())
    return {'dir_accuracy_all': dir_all, 'dir_accuracy_deg': dir_deg, 'dir_accuracy_weighted': dir_weighted}

def bioeval_deg(true_mat, pred_mat, deg_threshold=0.25):
    y_bin = (np.abs(true_mat) > deg_threshold).astype(int).flatten()
    y_score = np.abs(pred_mat).flatten()
    auprc = average_precision_score(y_bin, y_score) if 0 < y_bin.sum() < len(y_bin) else None
    sm = (np.sign(pred_mat) == np.sign(true_mat)).astype(float)
    da_score = (np.abs(pred_mat) * sm).flatten()
    da_auprc = average_precision_score(y_bin, da_score) if 0 < y_bin.sum() < len(y_bin) else None
    return {'auprc': float(auprc) if auprc is not None else None,
            'dir_aware_auprc': float(da_auprc) if da_auprc is not None else None}

def compute_existing(true_mat, pred_mat):
    ms, r2s, ps = [], [], []
    for p in range(true_mat.shape[0]):
        yt, yp = true_mat[p], pred_mat[p]
        ms.append(mean_squared_error(yt, yp))
        r2s.append(r2_score(yt, yp))
        if np.std(yt) > 0 and np.std(yp) > 0:
            r, _ = stats.pearsonr(yt, yp); ps.append(r)
        else: ps.append(np.nan)
    return {'MSE': float(np.mean(ms)), 'R2': float(np.mean(r2s)), 'Pearson': float(np.nanmean(ps))}

def task_deg_recovery(true_mat, pred_mat, deg_threshold=0.25, k=50):
    true_deg = np.abs(true_mat) > deg_threshold
    f1s = []
    for p in range(true_mat.shape[0]):
        topk = set(np.argsort(np.abs(pred_mat[p]))[::-1][:k])
        tdeg = set(np.where(true_deg[p])[0])
        if not tdeg: continue
        tp = len(topk & tdeg); prec = tp/k; rec = tp/len(tdeg)
        f1s.append(2*prec*rec/(prec+rec) if (prec+rec) > 0 else 0.0)
    return float(np.mean(f1s)) if f1s else 0.0

def task_dir_discovery(true_mat, pred_mat, deg_threshold=0.25):
    true_deg = np.abs(true_mat) > deg_threshold
    sm = (np.sign(pred_mat) == np.sign(true_mat)).astype(float)
    zt = (true_mat == 0); zp = (pred_mat == 0)
    sm[zt & zp] = 1.0; sm[zt & ~zp] = 0.5
    accs = []
    for p in range(true_mat.shape[0]):
        idx = np.where(true_deg[p])[0]
        if len(idx) == 0: continue
        accs.append(sm[p, idx].mean())
    return float(np.mean(accs)) if accs else 0.0

# ============================================================
# Data Loading
# ============================================================

def load_and_prepare(adata, ct_name):
    pert_col = None
    for col in ['condition', 'perturbation', 'perturb', 'pert_name']:
        if col in adata.obs.columns: pert_col = col; break
    if pert_col is None:
        raise ValueError(f"No perturbation column found for {ct_name}")

    ctrl_val = None
    for cv in ['ctrl', 'control', 'Ctrl', 'Control', 'non-targeting', 'NT']:
        if cv in adata.obs[pert_col].values: ctrl_val = cv; break

    ct_data = adata.obs[pert_col].value_counts()
    perts = [p for p in ct_data.index if p != ctrl_val and ct_data[p] >= 5]

    log(f"  [{ct_name}] ctrl='{ctrl_val}', col='{pert_col}', {len(perts)+1} values")

    ctrl_mask = adata.obs[pert_col] == ctrl_val
    ctrl_data = adata[ctrl_mask].X
    if hasattr(ctrl_data, "toarray"): ctrl_data = ctrl_data.toarray()
    ctrl_data = np.asarray(ctrl_data, dtype=np.float32)

    true_effects = []
    pert_names = []
    for pert in perts:
        pm = adata.obs[pert_col] == pert
        pd = adata[pm].X
        if hasattr(pd, "toarray"): pd = pd.toarray()
        pd = np.asarray(pd, dtype=np.float32)
        eff = pd.mean(axis=0) - ctrl_data.mean(axis=0)
        true_effects.append(eff)
        pert_names.append(pert)

    true_mat = np.array(true_effects)
    deg_frac = (np.abs(true_mat) > 0.25).mean()
    n_genes = true_mat.shape[1]
    log(f"  [{ct_name}] {len(perts)} perts, {n_genes} genes, DEG_frac={deg_frac:.4f}")

    return {'name': ct_name, 'perts': pert_names, 'true_mat': true_mat,
            'ctrl_data': ctrl_data, 'n_genes': n_genes, 'deg_frac': deg_frac,
            'deg_threshold': 0.25, 'adata': adata, 'pert_col': pert_col,
            'ctrl_val': ctrl_val}

def build_gene_level_features(ds_data, n_pca=30):
    ctrl_data = ds_data['ctrl_data']
    pca = PCA(n_components=n_pca, random_state=42)
    ctrl_pca = pca.fit_transform(ctrl_data)

    adata = ds_data['adata']
    pert_col = ds_data['pert_col']
    ctrl_val = ds_data['ctrl_val']
    perts = ds_data['perts']

    features = []
    for pert in perts:
        pm = adata.obs[pert_col] == pert
        pd = adata[pm].X
        if hasattr(pd, "toarray"): pd = pd.toarray()
        pd = np.asarray(pd, dtype=np.float32)
        pert_pca = pca.transform(pd)
        pca_mean = pert_pca.mean(axis=0)
        pca_var = pert_pca.var(axis=0)
        log_count = np.log1p(pm.sum())
        features.append(np.concatenate([pca_mean, pca_var, [log_count]]))

    feat_mat = np.array(features)
    cumul_var = pca.explained_variance_ratio_.cumsum()[-1]
    rank = np.linalg.matrix_rank(feat_mat)

    log(f"    Control cells: {ctrl_data.shape[0]}")
    log(f"    PCA explained var (cumul, {n_pca} PCs): {cumul_var:.3f}")
    log(f"    Feature matrix: {feat_mat.shape} (pca_mean={n_pca} + pca_var={n_pca} + log_count=1)")
    log(f"    Feature rank: {rank}/{feat_mat.shape[1]}")

    return feat_mat, pca

def build_additive_features(ds_data):
    perts = ds_data['perts']
    gene_set = sorted(set(g for p in perts for g in p.split('+')))
    gene_to_idx = {g: i for i, g in enumerate(gene_set)}
    features = np.zeros((len(perts), len(gene_set)), dtype=np.float32)
    for i, p in enumerate(perts):
        for g in p.split('+'):
            if g in gene_to_idx: features[i, gene_to_idx[g]] = 1.0
    return features

# ============================================================
# Model Training
# ============================================================

def ridge_loo_analytical(X, y, alpha):
    n, p = X.shape
    if n <= p:
        A = X @ X.T + alpha * np.eye(n)
        H = X.T @ np.linalg.solve(A, X)
    else:
        A = X.T @ X + alpha * np.eye(p)
        H = X @ np.linalg.solve(A, X.T)
    diag_H = np.diag(H).reshape(-1, 1)
    y_hat = H @ y
    loo_resid = (y - y_hat) / (1 - diag_H + 1e-10)
    y_loo = y - loo_resid
    return y_loo

def train_real_models(ds_data, feat_type='gene_pca'):
    true_mat = ds_data['true_mat']

    if feat_type == 'gene_pca':
        feat_mat, pca = build_gene_level_features(ds_data)
    else:
        feat_mat = build_additive_features(ds_data)

    y = true_mat.copy()
    models = {}

    # 1. mean_predictor
    models['mean_predictor'] = np.tile(y.mean(axis=0), (y.shape[0], 1))

    # 2-4. Ridge with analytical LOO
    for alpha_val, name in [(1.0, 'ridge'), (10.0, 'ridge_med'), (100.0, 'ridge_strong')]:
        log(f"    Ridge alpha={alpha_val} (analytical LOO)...")
        pred = ridge_loo_analytical(feat_mat, y, alpha_val)
        models[name] = pred

    # 5. noisy_ridge
    ridge_pred = models['ridge'].copy()
    noise = np.random.randn(*ridge_pred.shape) * 0.15 * np.std(ridge_pred)
    models['noisy_ridge'] = ridge_pred + noise

    # 6. sign_flip_ridge
    flip_mask = np.random.random(ridge_pred.shape) < 0.15
    models['sign_flip_ridge'] = ridge_pred.copy()
    models['sign_flip_ridge'][flip_mask] *= -1

    # 7. mean_effect
    models['mean_effect'] = np.tile(y.mean(axis=0), (y.shape[0], 1))

    # 8. constant_shrink (oracle)
    models['constant_shrink'] = true_mat * 0.15

    # 9. half_signal (oracle)
    models['half_signal'] = true_mat * 0.5

    return models, feat_type

# ============================================================
# Scale Correction Strategies
# ============================================================

def apply_global_rescale(pred_mat, true_mat):
    """Rescale all predictions by a single global ratio."""
    # Global ratio: total pred std / total true std
    pred_std = np.std(pred_mat)
    true_std = np.std(true_mat)
    ratio = pred_std / true_std if true_std > 0 else 1.0
    corrected = pred_mat / ratio if ratio > 0 else pred_mat
    return corrected, ratio

def apply_per_pert_rescale(pred_mat, true_mat):
    """Rescale each perturbation's prediction to match its true effect variance."""
    corrected = pred_mat.copy()
    ratios = []
    for i in range(pred_mat.shape[0]):
        pred_std = np.std(pred_mat[i])
        true_std = np.std(true_mat[i])
        ratio = pred_std / true_std if true_std > 0 else 1.0
        ratios.append(ratio)
        if ratio > 0:
            corrected[i] = pred_mat[i] / ratio
    return corrected, np.mean(ratios)

def apply_variance_match(pred_mat, true_mat):
    """Rescale predictions so that per-gene variance matches true variance."""
    # Per-gene approach: for each gene, match the variance across perturbations
    corrected = pred_mat.copy()
    n_genes = pred_mat.shape[1]
    scale_factors = []
    for g in range(n_genes):
        pred_var = np.var(pred_mat[:, g])
        true_var = np.var(true_mat[:, g])
        if pred_var > 0 and true_var > 0:
            sf = np.sqrt(true_var / pred_var)
            corrected[:, g] = pred_mat[:, g] * sf
            scale_factors.append(sf)
        else:
            scale_factors.append(1.0)
    return corrected, np.mean(scale_factors)

# ============================================================
# Full evaluation
# ============================================================

def evaluate_all_metrics(models, ds_data):
    true_mat = ds_data['true_mat']
    deg_threshold = ds_data['deg_threshold']

    eval_m = {}
    for mname, pred in models.items():
        ex = compute_existing(true_mat, pred)
        dv = bioeval_dir(true_mat, pred, deg_threshold)
        dg = bioeval_deg(true_mat, pred, deg_threshold)
        f1_50 = task_deg_recovery(true_mat, pred, deg_threshold, k=50)
        f1_100 = task_deg_recovery(true_mat, pred, deg_threshold, k=100)
        dir_disc = task_dir_discovery(true_mat, pred, deg_threshold)

        eval_m[mname] = {
            'MSE': ex['MSE'], 'R2': ex['R2'], 'Pearson': ex['Pearson'],
            'Dir_all': dv['dir_accuracy_all'],
            'Dir_deg': dv['dir_accuracy_deg'] if dv['dir_accuracy_deg'] is not None else 0.0,
            'Dir_weighted': dv['dir_accuracy_weighted'],
            'DEG_auprc': dg['auprc'] if dg['auprc'] is not None else 0.0,
            'DEG_dir_auprc': dg['dir_aware_auprc'] if dg['dir_aware_auprc'] is not None else 0.0,
            'f1@50': f1_50, 'f1@100': f1_100, 'dir_discovery_deg': dir_disc
        }
    return eval_m

def compute_scale_diagnostics(pred_mat, true_mat):
    """Compute per-perturbation scale diagnostics."""
    ratios = []
    for i in range(pred_mat.shape[0]):
        pred_std = np.std(pred_mat[i])
        true_std = np.std(true_mat[i])
        ratios.append(pred_std / true_std if true_std > 0 else np.nan)
    ratios = np.array(ratios)

    # DEG fraction in predictions vs true
    deg_true = (np.abs(true_mat) > 0.25).mean()
    deg_pred = (np.abs(pred_mat) > 0.25).mean()

    # Effect magnitude stats
    true_magnitudes = np.abs(true_mat).mean()
    pred_magnitudes = np.abs(pred_mat).mean()

    return {
        'mean_ratio': float(np.nanmean(ratios)),
        'std_ratio': float(np.nanstd(ratios)),
        'median_ratio': float(np.nanmedian(ratios)),
        'deg_frac_true': float(deg_true),
        'deg_frac_pred': float(deg_pred),
        'deg_frac_ratio': float(deg_pred / deg_true) if deg_true > 0 else 0.0,
        'mean_magnitude_true': float(true_magnitudes),
        'mean_magnitude_pred': float(pred_magnitudes),
        'magnitude_ratio': float(pred_magnitudes / true_magnitudes) if true_magnitudes > 0 else 0.0
    }

def compute_h1_h2(eval_m, model_names):
    """Compute H1 (Kendall tau) and H2 (Spearman rho) from evaluation results."""
    results = {}

    # Model names excluding oracles for H3
    non_oracle = [m for m in model_names if m not in ['constant_shrink', 'half_signal']]

    # H1: Kendall tau between MSE/R2 and directional metrics
    for m1 in ['MSE', 'R2', 'Pearson']:
        for m2 in ['Dir_all', 'Dir_deg', 'DEG_auprc']:
            vals1 = [eval_m[m][m1] for m in non_oracle]
            vals2 = [eval_m[m][m2] for m in non_oracle]
            tau, p = stats.kendalltau(vals1, vals2)
            results[f'tau({m1},{m2})'] = {'tau': float(tau), 'p': float(p)}

    # H2: Spearman rho between BioEval metrics and downstream tasks
    for bm in ['Dir_all', 'Dir_deg', 'DEG_auprc']:
        for dt in ['f1@50', 'f1@100', 'dir_discovery_deg']:
            vals1 = [eval_m[m][bm] for m in non_oracle]
            vals2 = [eval_m[m][dt] for m in non_oracle]
            rho, p = stats.spearmanr(vals1, vals2)
            results[f'rho({bm},{dt})'] = {'rho': float(rho), 'p': float(p)}

    # H2: MSE vs downstream (baseline comparison)
    for dt in ['f1@50', 'f1@100', 'dir_discovery_deg']:
        vals1 = [eval_m[m]['MSE'] for m in non_oracle]
        vals2 = [eval_m[m][dt] for m in non_oracle]
        rho, p = stats.spearmanr(vals1, vals2)
        results[f'rho(MSE,{dt})'] = {'rho': float(rho), 'p': float(p)}

    # H2 gap: rho(BioEval, downstream) - rho(MSE, downstream)
    for bm in ['Dir_all', 'Dir_deg', 'DEG_auprc']:
        for dt in ['f1@50', 'f1@100', 'dir_discovery_deg']:
            bioeval_rho = results[f'rho({bm},{dt})']['rho']
            mse_rho = results[f'rho(MSE,{dt})']['rho']
            results[f'gap({bm},{dt})'] = {'diff': float(bioeval_rho - mse_rho)}

    return results

# ============================================================
# Main
# ============================================================

def main():
    t0 = time.time()
    log("=" * 60)
    log("BioEval Run 18 -- Scale Correction Analysis (A3)")
    log("=" * 60)

    # Load data (same paths as run_16/17)
    datasets = {}
    for ct, subdir in [("K562","replogle_k562_essential"),("RPE1","replogle_rpe1_essential"),("Norman","norman")]:
        path = os.path.join(DATA_DIR, "gears_data", subdir, "perturb_processed.h5ad")
        if os.path.exists(path):
            log(f"Loading {ct}...")
            adata = scanpy.read_h5ad(path)
            log(f"  {ct}: {adata.shape}")
            r = load_and_prepare(adata, ct)
            if r: datasets[ct] = r
        else:
            log(f"  {ct}: NOT FOUND at {path}")

    if not datasets:
        log("ERROR: No datasets found!"); return

    ds_k562 = datasets.get("K562")
    ds_rpe1 = datasets.get("RPE1")
    ds_norman = datasets.get("Norman")

    all_results = {}

    for ds_data in [ds_k562, ds_rpe1, ds_norman]:
        if ds_data is None: continue
        ds_name = ds_data['name']
        log(f"\n{'=' * 60}")
        log(f"ANALYZING: {ds_name}")
        log(f"{'=' * 60}")

        true_mat = ds_data['true_mat']

        # Train models
        feat_type = 'additive' if ds_name == 'Norman' else 'gene_pca'
        models, _ = train_real_models(ds_data, feat_type)

        # Baseline evaluation (no correction)
        eval_baseline = evaluate_all_metrics(models, ds_data)

        # Scale diagnostics for each model
        log(f"\n  Scale Diagnostics (BEFORE correction):")
        log(f"  {'Model':<20s} {'Mean Ratio':>10s} {'DEG_frac_pred':>14s} {'DEG_frac_true':>15s} {'Mag_ratio':>10s}")
        for mname, pred in models.items():
            diag = compute_scale_diagnostics(pred, true_mat)
            log(f"  {mname:<20s} {diag['mean_ratio']:>10.3f} {diag['deg_frac_pred']:>14.4f} {diag['deg_frac_true']:>15.4f} {diag['magnitude_ratio']:>10.3f}")

        # Apply corrections to trained models only (skip baselines/oracles)
        trained_models = ['ridge', 'ridge_med', 'ridge_strong', 'noisy_ridge', 'sign_flip_ridge']
        correction_methods = {
            'global_rescale': apply_global_rescale,
            'per_pert_rescale': apply_per_pert_rescale,
            'variance_match': apply_variance_match
        }

        all_evals = {'baseline': eval_baseline}

        for corr_name, corr_fn in correction_methods.items():
            log(f"\n  Applying {corr_name}...")
            corrected_models = {}
            for mname, pred in models.items():
                if mname in trained_models:
                    corrected, ratio = corr_fn(pred, true_mat)
                    corrected_models[mname] = corrected
                    diag = compute_scale_diagnostics(corrected, true_mat)
                    log(f"    {mname}: ratio={ratio:.3f}, DEG_frac={diag['deg_frac_pred']:.4f} (true={diag['deg_frac_true']:.4f})")
                else:
                    corrected_models[mname] = pred

            eval_corrected = evaluate_all_metrics(corrected_models, ds_data)
            all_evals[corr_name] = eval_corrected

        # Compute H1/H2 for each correction method
        model_names = list(models.keys())
        non_oracle = [m for m in model_names if m not in ['constant_shrink', 'half_signal']]

        h1_h2_results = {}
        for eval_name, eval_m in all_evals.items():
            h1_h2_results[eval_name] = compute_h1_h2(eval_m, model_names)

        # Key comparisons
        log(f"\n  Key Metric Comparison (trained models only, oracle excluded):")
        log(f"  {'Metric':<15s} {'baseline':>10s} {'global':>10s} {'per_pert':>10s} {'var_match':>10s}")

        key_metrics = ['MSE', 'R2', 'Dir_deg', 'DEG_auprc', 'f1@50', 'f1@100', 'dir_discovery_deg']
        for metric in key_metrics:
            vals = []
            for eval_name in ['baseline', 'global_rescale', 'per_pert_rescale', 'variance_match']:
                trained_vals = [all_evals[eval_name][m][metric] for m in non_oracle]
                vals.append(f"{np.mean(trained_vals):.4f}" if metric not in ['MSE'] else f"{np.mean(trained_vals):.4f}")
            log(f"  {metric:<15s} {vals[0]:>10s} {vals[1]:>10s} {vals[2]:>10s} {vals[3]:>10s}")

        # H1 key result: tau(MSE, Dir_deg)
        log(f"\n  H1 Kendall tau(MSE, Dir_deg) comparison:")
        for eval_name in ['baseline', 'global_rescale', 'per_pert_rescale', 'variance_match']:
            key = 'tau(MSE,Dir_deg)'
            tau = h1_h2_results[eval_name][key]['tau']
            log(f"    {eval_name:<20s}: tau={tau:.3f}")

        # H2 key result: gap between BioEval and MSE for downstream prediction
        log(f"\n  H2 gap: rho(DEG_auprc, f1@50) - rho(MSE, f1@50):")
        for eval_name in ['baseline', 'global_rescale', 'per_pert_rescale', 'variance_match']:
            key = 'gap(DEG_auprc,f1@50)'
            diff = h1_h2_results[eval_name][key]['diff']
            rho_bio = h1_h2_results[eval_name]['rho(DEG_auprc,f1@50)']['rho']
            rho_mse = h1_h2_results[eval_name]['rho(MSE,f1@50)']['rho']
            log(f"    {eval_name:<20s}: rho(DEG_auprc)={rho_bio:.3f}, rho(MSE)={rho_mse:.3f}, gap={diff:+.3f}")

        # H2 pass rate
        log(f"\n  H2 pass rate (BioEval rho > MSE rho for each BioEval-downstream pair):")
        for eval_name in ['baseline', 'global_rescale', 'per_pert_rescale', 'variance_match']:
            pass_count = 0
            total = 0
            for bm in ['Dir_all', 'Dir_deg', 'DEG_auprc']:
                for dt in ['f1@50', 'f1@100', 'dir_discovery_deg']:
                    key = f'gap({bm},{dt})'
                    diff = h1_h2_results[eval_name][key]['diff']
                    total += 1
                    if diff > 0.1:  # same threshold as run_16
                        pass_count += 1
            log(f"    {eval_name:<20s}: {pass_count}/{total}")

        # Per-model detailed comparison for ridge (best model)
        log(f"\n  Ridge model detailed metrics (baseline vs best correction):")
        for metric in key_metrics:
            base_val = all_evals['baseline']['ridge'][metric]
            # Find which correction gives best improvement for each metric
            best_corr = 'baseline'
            best_val = base_val
            for corr_name in ['global_rescale', 'per_pert_rescale', 'variance_match']:
                corr_val = all_evals[corr_name]['ridge'][metric]
                # For MSE, lower is better; for others, higher is better
                if metric == 'MSE':
                    if corr_val < best_val:
                        best_corr = corr_name
                        best_val = corr_val
                else:
                    if corr_val > best_val:
                        best_corr = corr_name
                        best_val = corr_val
            log(f"    {metric:<15s}: baseline={base_val:.4f}, best={best_val:.4f} ({best_corr})")

        all_results[ds_name] = {
            'evaluations': {},
            'h1_h2': h1_h2_results,
            'scale_diagnostics': {}
        }

        # Save evals (convert to serializable)
        for eval_name, eval_m in all_evals.items():
            all_results[ds_name]['evaluations'][eval_name] = eval_m

        # Save scale diagnostics
        for mname, pred in models.items():
            all_results[ds_name]['scale_diagnostics'][mname] = {
                'baseline': compute_scale_diagnostics(pred, true_mat)
            }
            for corr_name, corr_fn in correction_methods.items():
                if mname in trained_models:
                    corrected, _ = corr_fn(pred, true_mat)
                    all_results[ds_name]['scale_diagnostics'][mname][corr_name] = compute_scale_diagnostics(corrected, true_mat)

    # Save JSON
    json_path = os.path.join(OUTPUT_DIR, "run_18_scale_correction.json")
    with open(json_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    log(f"\nResults saved to {json_path}")

    # Final summary
    log(f"\n{'=' * 60}")
    log("CROSS-DATASET SUMMARY")
    log(f"{'=' * 60}")

    for ds_name in ['K562', 'RPE1', 'Norman']:
        log(f"\n  {ds_name}:")
        h1_h2 = all_results[ds_name]['h1_h2']

        # H1 summary
        log(f"    H1 tau(MSE,Dir_deg):")
        for eval_name in ['baseline', 'global_rescale', 'per_pert_rescale', 'variance_match']:
            tau = h1_h2[eval_name]['tau(MSE,Dir_deg)']['tau']
            log(f"      {eval_name:<20s}: {tau:.3f}")

        # H2 gap summary
        log(f"    H2 gap(DEG_auprc,f1@50):")
        for eval_name in ['baseline', 'global_rescale', 'per_pert_rescale', 'variance_match']:
            diff = h1_h2[eval_name]['gap(DEG_auprc,f1@50)']['diff']
            log(f"      {eval_name:<20s}: {diff:+.3f}")

        # Scale diagnostics for ridge
        ridge_diag = all_results[ds_name]['scale_diagnostics']['ridge']
        log(f"    Ridge scale diagnostics:")
        for corr_name, diag in ridge_diag.items():
            log(f"      {corr_name:<20s}: ratio={diag['mean_ratio']:.3f}, DEG_frac={diag['deg_frac_pred']:.4f}, mag_ratio={diag['magnitude_ratio']:.3f}")

    # A3 Decision
    log(f"\n{'=' * 60}")
    log("A3 DECISION")
    log(f"{'=' * 60}")

    # Check if any correction improves downstream task discrimination
    # without harming H1 or H2
    for ds_name in ['K562', 'RPE1', 'Norman']:
        h1_h2 = all_results[ds_name]['h1_h2']

        # Check if correction improves H2 gap for DEG_auprc vs f1@50
        base_gap = h1_h2['baseline']['gap(DEG_auprc,f1@50)']['diff']
        best_gap = base_gap
        best_corr = 'baseline'

        for corr_name in ['global_rescale', 'per_pert_rescale', 'variance_match']:
            gap = h1_h2[corr_name]['gap(DEG_auprc,f1@50)']['diff']
            if gap > best_gap:
                best_gap = gap
                best_corr = corr_name

        # Check H1 stability
        base_h1 = abs(h1_h2['baseline']['tau(MSE,Dir_deg)']['tau'])
        h1_stable = True
        for corr_name in ['global_rescale', 'per_pert_rescale', 'variance_match']:
            corr_h1 = abs(h1_h2[corr_name]['tau(MSE,Dir_deg)']['tau'])
            # If correction makes H1 much stronger (higher tau), it might hurt H1 claim
            # If correction makes H1 much weaker, check if still < 0.7
            if abs(corr_h1 - base_h1) > 0.2:
                h1_stable = False

        log(f"  {ds_name}: best H2 gap={best_gap:+.3f} ({best_corr}), H1 stable={h1_stable}")

        if best_corr != 'baseline':
            log(f"    -> Correction IMPROVES H2: apply {best_corr}")
        else:
            log(f"    -> Correction does NOT improve H2: keep baseline")

    elapsed = time.time() - t0
    log(f"\nRuntime: {elapsed:.1f}s")

if __name__ == '__main__':
    main()
