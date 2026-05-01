"""
BioEval Run 17 -- Bootstrap CI for all tau/rho estimates
Computes bootstrap confidence intervals for H1 (Kendall tau) and H2 (Spearman rho)
across all 3 datasets using the gene PCA Ridge results from run_16.

Strategy:
- Load run_16 model results (9 models per dataset)
- Bootstrap resample models (with replacement) B=10000 times
- For each bootstrap sample, recompute tau and rho
- Report 95% CI (2.5th and 97.5th percentiles)
- Also compute p-values via bootstrap (fraction of samples with opposite sign)
"""
import json, os, warnings, time
import numpy as np
import scanpy
from scipy import stats
from sklearn.metrics import mean_squared_error, r2_score, average_precision_score
from sklearn.linear_model import Ridge
from sklearn.decomposition import PCA

warnings.filterwarnings('ignore')

OUTPUT_DIR = r"C:\test\llm_project\nexus-science-win\outputs\analysis\run_17"
DATA_DIR = r"C:\test\llm_project\nexus-science-win\outputs\analysis\run_04\data"
os.makedirs(OUTPUT_DIR, exist_ok=True)

B_BOOT = 10000  # number of bootstrap samples
np.random.seed(42)

def log(msg):
    try:
        print(msg, flush=True)
    except UnicodeEncodeError:
        print(msg.encode('ascii', errors='replace').decode('ascii'), flush=True)

# ============================================================
# BioEval Metrics (same as run_16)
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
# Data Loading (same as run_16)
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
# Model Training (same as run_16)
# ============================================================

def ridge_loo_analytical(X, y, alpha):
    n, p = X.shape
    if n <= p:
        A = X @ X.T + alpha * np.eye(n)
        H = X.T @ np.linalg.solve(A, X)
    else:
        A = X.T @ X + alpha * np.eye(p)
        H = X @ np.linalg.solve(A, X.T)
    diag_H = np.diag(H).reshape(-1, 1)  # (n, 1) for broadcasting with (n, p)
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
# Full evaluation (returns dict of metric values per model)
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

# ============================================================
# Bootstrap CI for H1 (Kendall tau) and H2 (Spearman rho)
# ============================================================

def bootstrap_h1_tau(eval_m, metric_pairs, B=B_BOOT):
    """
    Bootstrap CI for Kendall tau between metric pairs.
    Resample models (N=9) with replacement, recompute tau.
    """
    model_names = [m for m in eval_m.keys()]
    results = {}

    for m1_name, m2_name in metric_pairs:
        m1_vals = np.array([eval_m[m][m1_name] for m in model_names])
        m2_vals = np.array([eval_m[m][m2_name] for m in model_names])

        # Observed tau
        tau_obs, p_obs = stats.kendalltau(m1_vals, m2_vals)

        # Bootstrap
        tau_boot = np.zeros(B)
        for b in range(B):
            idx = np.random.choice(len(model_names), size=len(model_names), replace=True)
            tau_b, _ = stats.kendalltau(m1_vals[idx], m2_vals[idx])
            tau_boot[b] = tau_b if not np.isnan(tau_b) else 0.0

        ci_lo = float(np.percentile(tau_boot, 2.5))
        ci_hi = float(np.percentile(tau_boot, 97.5))
        se = float(np.std(tau_boot))

        # Bootstrap p-value: fraction where tau has opposite sign to observed
        if tau_obs > 0:
            p_boot = float((tau_boot <= 0).mean())
        elif tau_obs < 0:
            p_boot = float((tau_boot >= 0).mean())
        else:
            p_boot = 1.0

        results[f"tau({m1_name},{m2_name})"] = {
            'observed': float(tau_obs),
            'CI_95': [ci_lo, ci_hi],
            'SE': se,
            'p_asymptotic': float(p_obs),
            'p_bootstrap': p_boot,
            'significant': ci_lo > 0 or ci_hi < 0  # CI doesn't include 0
        }

    return results

def bootstrap_h2_rho(eval_m, metric_pairs, B=B_BOOT):
    """
    Bootstrap CI for Spearman rho between metric values and downstream task values.
    Resample models with replacement, recompute rho.
    """
    model_names = [m for m in eval_m.keys()]
    results = {}

    for metric_name, ds_name in metric_pairs:
        metric_vals = np.array([eval_m[m][metric_name] for m in model_names])
        ds_vals = np.array([eval_m[m][ds_name] for m in model_names])

        # Observed rho
        rho_obs, p_obs = stats.spearmanr(metric_vals, ds_vals)

        # Bootstrap
        rho_boot = np.zeros(B)
        for b in range(B):
            idx = np.random.choice(len(model_names), size=len(model_names), replace=True)
            rho_b, _ = stats.spearmanr(metric_vals[idx], ds_vals[idx])
            rho_boot[b] = rho_b if not np.isnan(rho_b) else 0.0

        ci_lo = float(np.percentile(rho_boot, 2.5))
        ci_hi = float(np.percentile(rho_boot, 97.5))
        se = float(np.std(rho_boot))

        if rho_obs > 0:
            p_boot = float((rho_boot <= 0).mean())
        elif rho_obs < 0:
            p_boot = float((rho_boot >= 0).mean())
        else:
            p_boot = 1.0

        results[f"rho({metric_name},{ds_name})"] = {
            'observed': float(rho_obs),
            'CI_95': [ci_lo, ci_hi],
            'SE': se,
            'p_asymptotic': float(p_obs),
            'p_bootstrap': p_boot,
            'significant': ci_lo > 0 or ci_hi < 0
        }

    return results

def bootstrap_h2_diff(eval_m, diff_pairs, B=B_BOOT):
    """
    Bootstrap CI for rho difference: rho(bioeval_metric, ds) - rho(MSE, ds).
    This is the key H2 test: is BioEval significantly better than MSE?
    """
    model_names = [m for m in eval_m.keys()]
    results = {}

    for bioeval_metric, ds_name in diff_pairs:
        be_vals = np.array([eval_m[m][bioeval_metric] for m in model_names])
        mse_vals = np.array([eval_m[m]['MSE'] for m in model_names])
        ds_vals = np.array([eval_m[m][ds_name] for m in model_names])

        rho_be_obs, _ = stats.spearmanr(be_vals, ds_vals)
        rho_mse_obs, _ = stats.spearmanr(mse_vals, ds_vals)
        diff_obs = rho_be_obs - rho_mse_obs

        # Bootstrap
        diff_boot = np.zeros(B)
        for b in range(B):
            idx = np.random.choice(len(model_names), size=len(model_names), replace=True)
            rho_be_b, _ = stats.spearmanr(be_vals[idx], ds_vals[idx])
            rho_mse_b, _ = stats.spearmanr(mse_vals[idx], ds_vals[idx])
            diff_boot[b] = (rho_be_b - rho_mse_b) if not (np.isnan(rho_be_b) or np.isnan(rho_mse_b)) else 0.0

        ci_lo = float(np.percentile(diff_boot, 2.5))
        ci_hi = float(np.percentile(diff_boot, 97.5))
        se = float(np.std(diff_boot))
        p_boot = float((diff_boot <= 0).mean())  # H0: diff <= 0

        results[f"diff_rho({bioeval_metric},{ds_name})-rho(MSE,{ds_name})"] = {
            'rho_bioeval_obs': float(rho_be_obs),
            'rho_mse_obs': float(rho_mse_obs),
            'diff_observed': float(diff_obs),
            'CI_95': [ci_lo, ci_hi],
            'SE': se,
            'p_bootstrap': p_boot,
            'significant': ci_lo > 0  # CI doesn't include 0 (one-sided: bioeval > MSE)
        }

    return results

# ============================================================
# Main
# ============================================================

def main():
    t0 = time.time()
    log("=" * 60)
    log("BioEval Run 17 -- Bootstrap CI for tau/rho")
    log("=" * 60)

    # Load data
    log("Loading K562...")
    adata_k562 = scanpy.read_h5ad(os.path.join(DATA_DIR, "gears_data", "replogle_k562_essential", "perturb_processed.h5ad"))
    ds_k562 = load_and_prepare(adata_k562, "K562")

    log("Loading RPE1...")
    adata_rpe1 = scanpy.read_h5ad(os.path.join(DATA_DIR, "gears_data", "replogle_rpe1_essential", "perturb_processed.h5ad"))
    ds_rpe1 = load_and_prepare(adata_rpe1, "RPE1")

    log("Loading Norman...")
    adata_norman = scanpy.read_h5ad(os.path.join(DATA_DIR, "gears_data", "norman", "perturb_processed.h5ad"))
    ds_norman = load_and_prepare(adata_norman, "Norman")

    datasets = [
        (ds_k562, 'gene_pca'),
        (ds_rpe1, 'gene_pca'),
        (ds_norman, 'additive'),
    ]

    all_results = {}

    for ds_data, feat_type in datasets:
        ds_name = ds_data['name']
        log("")
        log("=" * 60)
        log(f"ANALYZING: {ds_name} -- Bootstrap CI")
        log("=" * 60)

        # Train models (same as run_16)
        models, _ = train_real_models(ds_data, feat_type)
        eval_m = evaluate_all_metrics(models, ds_data)

        model_names = list(eval_m.keys())
        log(f"  Models ({len(model_names)}): {model_names}")

        # ---- H1: Bootstrap CI for Kendall tau ----
        log("\n  H1: Bootstrap CI for Kendall tau...")
        h1_pairs = [
            ('MSE', 'Dir_all'),
            ('MSE', 'Dir_deg'),
            ('MSE', 'DEG_auprc'),
            ('R2', 'Dir_deg'),
            ('Pearson', 'Dir_deg'),
        ]
        h1_results = bootstrap_h1_tau(eval_m, h1_pairs)

        for key, val in h1_results.items():
            sig_mark = "***" if val['significant'] else ""
            log(f"    {key}: tau={val['observed']:.3f} CI=[{val['CI_95'][0]:.3f}, {val['CI_95'][1]:.3f}] "
                f"SE={val['SE']:.3f} p_boot={val['p_bootstrap']:.4f} {sig_mark}")

        # ---- H2: Bootstrap CI for Spearman rho ----
        log("\n  H2: Bootstrap CI for Spearman rho...")
        h2_pairs = [
            ('Dir_all', 'f1@50'),
            ('Dir_deg', 'f1@50'),
            ('Dir_weighted', 'f1@50'),
            ('DEG_auprc', 'f1@50'),
            ('DEG_dir_auprc', 'f1@50'),
            ('Dir_all', 'f1@100'),
            ('Dir_deg', 'f1@100'),
            ('DEG_auprc', 'f1@100'),
            ('Dir_all', 'dir_discovery_deg'),
            ('Dir_deg', 'dir_discovery_deg'),
            ('DEG_auprc', 'dir_discovery_deg'),
        ]
        h2_results = bootstrap_h2_rho(eval_m, h2_pairs)

        for key, val in h2_results.items():
            sig_mark = "***" if val['significant'] else ""
            log(f"    {key}: rho={val['observed']:.3f} CI=[{val['CI_95'][0]:.3f}, {val['CI_95'][1]:.3f}] "
                f"SE={val['SE']:.3f} p_boot={val['p_bootstrap']:.4f} {sig_mark}")

        # ---- H2: Bootstrap CI for rho DIFF (BioEval - MSE) ----
        log("\n  H2: Bootstrap CI for rho DIFF (BioEval - MSE)...")
        diff_pairs = [
            ('Dir_all', 'f1@50'),
            ('Dir_deg', 'f1@50'),
            ('DEG_auprc', 'f1@50'),
            ('Dir_all', 'f1@100'),
            ('Dir_deg', 'f1@100'),
            ('DEG_auprc', 'f1@100'),
            ('Dir_all', 'dir_discovery_deg'),
            ('Dir_deg', 'dir_discovery_deg'),
            ('DEG_auprc', 'dir_discovery_deg'),
        ]
        diff_results = bootstrap_h2_diff(eval_m, diff_pairs)

        for key, val in diff_results.items():
            sig_mark = "***" if val['significant'] else ""
            log(f"    {key}: diff={val['diff_observed']:.3f} CI=[{val['CI_95'][0]:.3f}, {val['CI_95'][1]:.3f}] "
                f"p_boot={val['p_bootstrap']:.4f} {sig_mark}")

        all_results[ds_name] = {
            'h1_tau': h1_results,
            'h2_rho': h2_results,
            'h2_diff': diff_results,
        }

    # ---- Cross-dataset summary ----
    log("\n" + "=" * 60)
    log("CROSS-DATASET SUMMARY")
    log("=" * 60)

    for ds_name in ['K562', 'RPE1', 'Norman']:
        r = all_results[ds_name]
        log(f"\n  {ds_name}:")

        # H1 summary
        log("  H1 (Kendall tau):")
        reversal_count = 0
        for key, val in r['h1_tau'].items():
            ci = val['CI_95']
            if ci[0] < 0 and ci[1] > 0:
                reversal_count += 1
            sig = "SIG" if val['significant'] else "ns"
            log(f"    {key}: {val['observed']:.3f} [{ci[0]:.3f}, {ci[1]:.3f}] {sig}")
        log(f"    CI includes 0 (reversal possible): {reversal_count}/{len(r['h1_tau'])}")

        # H2 diff summary
        log("  H2 (rho diff BioEval - MSE):")
        sig_count = 0
        for key, val in r['h2_diff'].items():
            if val['significant']: sig_count += 1
            ci = val['CI_95']
            sig = "SIG" if val['significant'] else "ns"
            log(f"    {key}: diff={val['diff_observed']:.3f} [{ci[0]:.3f}, {ci[1]:.3f}] {sig}")
        log(f"    Significant (CI excludes 0): {sig_count}/{len(r['h2_diff'])}")

    # Save results
    # Convert numpy types for JSON serialization
    def convert(obj):
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, dict): return {k: convert(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)): return [convert(x) for x in obj]
        if isinstance(obj, bool): return obj
        return obj

    save_path = os.path.join(OUTPUT_DIR, "run_17_bootstrap_ci.json")
    with open(save_path, 'w') as f:
        json.dump(convert(all_results), f, indent=2)
    log(f"\nResults saved to {save_path}")

    elapsed = time.time() - t0
    log(f"Runtime: {elapsed:.1f}s")

if __name__ == "__main__":
    main()
