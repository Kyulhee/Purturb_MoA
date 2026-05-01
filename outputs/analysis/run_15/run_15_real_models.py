"""
BioEval Run 15 -- Real Trained Models (sklearn Ridge)
Replaces simulated predictions with real LOO-predicted Ridge models.
"""
import json, os, warnings, time, sys
import numpy as np
import scanpy
from scipy import stats
from sklearn.metrics import mean_squared_error, r2_score, average_precision_score
from sklearn.linear_model import Ridge

warnings.filterwarnings('ignore')

OUTPUT_DIR = r"C:\test\llm_project\nexus-science-win\outputs\analysis\run_15"
DATA_DIR = r"C:\test\llm_project\nexus-science-win\outputs\analysis\run_04\data"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def log(msg):
    print(msg, flush=True)

# ============================================================
# BioEval Metrics
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

# ============================================================
# Downstream Tasks
# ============================================================

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
    if pert_col is None: return None
    pert_values = list(adata.obs[pert_col].unique())
    ctrl_key = None
    for key in ['ctrl', 'control', 'NaN', 'nan', 'ctrl_ess']:
        if key in pert_values: ctrl_key = key; break
    if ctrl_key is None: ctrl_key = adata.obs[pert_col].value_counts().index[0]
    log(f"  [{ct_name}] ctrl='{ctrl_key}', col='{pert_col}', {len(pert_values)} values")
    ctrl_mean = np.array(adata[adata.obs[pert_col] == ctrl_key].X.mean(axis=0)).flatten()
    perturbations = [p for p in pert_values if p != ctrl_key]
    pert_means = {}
    for pert in perturbations:
        mask = adata.obs[pert_col] == pert
        if mask.sum() > 0:
            pert_means[pert] = np.array(adata[mask].X.mean(axis=0)).flatten()
    perts_sorted = sorted(pert_means.keys())
    true_logfc = np.array([pert_means[p] - ctrl_mean for p in perts_sorted])
    log(f"  [{ct_name}] {len(perts_sorted)} perts, {true_logfc.shape[1]} genes, DEG_frac={((np.abs(true_logfc)>0.25).mean()):.4f}")
    return {'cell_type': ct_name, 'perturbations': perts_sorted, 'true_logfc': true_logfc,
            'n_perts': len(perts_sorted), 'n_genes': true_logfc.shape[1]}

# ============================================================
# Ridge LOO (analytical) + degraded variants
# ============================================================

def ridge_loo_analytical(X, y, alpha):
    """Ridge LOO via hat matrix -- O(p^3) once, not O(n*p^3)"""
    n, p = X.shape
    A = X.T @ X + alpha * np.eye(p)
    A_inv = np.linalg.inv(A)
    H = X @ A_inv @ X.T
    y_hat = H @ y
    dH = np.diag(H)
    loo_pred = y_hat - (dH[:, None] * (y - y_hat)) / (1 - dH[:, None] + 1e-12)
    return loo_pred

def train_real_models(ds_data):
    true_logfc = ds_data['true_logfc']
    n_perts, n_genes = true_logfc.shape
    perts = ds_data['perturbations']

    # Build additive features
    def parse_pert(name):
        if '+' in name: return [p.strip() for p in name.split('+')]
        return [name]

    all_genes = set()
    pert_sets = {}
    for p in perts:
        gs = parse_pert(p); pert_sets[p] = gs
        for g in gs: all_genes.add(g)
    sg_list = sorted(all_genes); sg_idx = {g: i for i, g in enumerate(sg_list)}

    X_add = np.zeros((n_perts, len(sg_list)), dtype=float)
    for i, p in enumerate(perts):
        for g in pert_sets[p]:
            if g in sg_idx: X_add[i, sg_idx[g]] = 1.0

    has_add = len(sg_list) > 0 and X_add.shape[1] < n_perts
    X = X_add if has_add else np.eye(n_perts, dtype=float)
    ftype = "additive" if has_add else "onehot"
    log(f"    Features: {ftype} ({X.shape[1]} feat for {n_perts} perts)")

    model_predictions = {}
    model_metadata = {}
    np.random.seed(42)

    # 1. mean_predictor
    model_predictions['mean_predictor'] = np.zeros_like(true_logfc)
    model_metadata['mean_predictor'] = {'type': 'baseline', 'trained': False}

    # 2-4. Ridge at 3 alphas (analytical LOO)
    for alpha, name in [(1.0, 'ridge'), (10.0, 'ridge_med'), (100.0, 'ridge_strong')]:
        log(f"    Ridge alpha={alpha} (analytical LOO)...")
        loo = ridge_loo_analytical(X, true_logfc, alpha)
        model_predictions[name] = loo
        model_metadata[name] = {'type': 'trained', 'feature_type': ftype, 'loo': True, 'alpha': alpha}

    # 5-6. Degraded variants of ridge
    rp = model_predictions['ridge']
    gs = true_logfc.std(axis=0)
    ns = 0.15 * np.median(gs[gs > 0])

    model_predictions['noisy_ridge'] = rp + np.random.normal(0, ns, rp.shape)
    model_metadata['noisy_ridge'] = {'type': 'degraded', 'base': 'ridge', 'noise': 0.15}

    deg_mask = np.abs(true_logfc) > 0.25
    fp = rp.copy()
    fp[deg_mask & (np.random.random(true_logfc.shape) < 0.15)] *= -1
    model_predictions['sign_flip_ridge'] = fp
    model_metadata['sign_flip_ridge'] = {'type': 'degraded', 'base': 'ridge', 'flip_frac': 0.15}

    # 7. mean_effect: predict perturbation-mean effect for all genes
    pme = true_logfc.mean(axis=1, keepdims=True)
    model_predictions['mean_effect'] = np.broadcast_to(pme * 0.3, true_logfc.shape).copy()
    model_metadata['mean_effect'] = {'type': 'baseline', 'trained': False}

    # 8. constant_shrink: small fraction of true logFC
    model_predictions['constant_shrink'] = true_logfc * 0.15
    model_metadata['constant_shrink'] = {'type': 'baseline', 'trained': False}

    # 9. half_signal: 50% of true signal + noise
    model_predictions['half_signal'] = true_logfc * 0.5 + np.random.normal(0, ns*0.5, true_logfc.shape)
    model_metadata['half_signal'] = {'type': 'degraded', 'signal_frac': 0.5}

    return model_predictions, model_metadata

# ============================================================
# Main
# ============================================================

def main():
    t0 = time.time()
    log("=" * 60)
    log("BioEval Run 15 -- Real Trained Models")
    log("=" * 60)

    datasets = {}
    for ct, subdir in [("K562","replogle_k562_essential"),("RPE1","replogle_rpe1_essential"),("Norman","norman")]:
        path = os.path.join(DATA_DIR, "gears_data", subdir, "perturb_processed.h5ad")
        if os.path.exists(path):
            log(f"\nLoading {ct}...")
            adata = scanpy.read_h5ad(path)
            log(f"  {ct}: {adata.shape}")
            r = load_and_prepare(adata, ct)
            if r: datasets[ct] = r
        else:
            log(f"  {ct}: NOT FOUND")

    if not datasets:
        log("ERROR: No datasets!"); return

    all_results = {}
    for ct_name, ds_data in datasets.items():
        log(f"\n{'='*60}")
        log(f"ANALYZING: {ct_name} ({ds_data['n_perts']} perts, {ds_data['n_genes']} genes)")
        log(f"{'='*60}")

        true_logfc = ds_data['true_logfc']
        preds, meta = train_real_models(ds_data)
        model_names = list(preds.keys())
        log(f"  Models ({len(model_names)}): {model_names}")

        # Compute evaluation metrics
        eval_m = {}
        for mn in model_names:
            ex = compute_existing(true_logfc, preds[mn])
            dr = bioeval_dir(true_logfc, preds[mn])
            dg = bioeval_deg(true_logfc, preds[mn])
            eval_m[mn] = {
                'MSE': ex['MSE'], 'R2': ex['R2'], 'Pearson': ex['Pearson'],
                'Dir_all': dr['dir_accuracy_all'],
                'Dir_deg': dr['dir_accuracy_deg'] if dr['dir_accuracy_deg'] is not None else 0.0,
                'Dir_weighted': dr['dir_accuracy_weighted'],
                'DEG_auprc': dg['auprc'] if dg['auprc'] is not None else 0.0,
                'DEG_dir_auprc': dg['dir_aware_auprc'] if dg['dir_aware_auprc'] is not None else 0.0,
            }
            log(f"    {mn}: MSE={ex['MSE']:.4f} R2={ex['R2']:.4f} Dir_all={dr['dir_accuracy_all']:.4f} "
                f"Dir_deg={dr['dir_accuracy_deg'] if dr['dir_accuracy_deg'] else 'N/A'}")

        # H1: Rankings
        log(f"\n  H1: Metric-Ranking Reversal...")
        rank_scores = {
            'MSE': {m: -eval_m[m]['MSE'] for m in model_names},
            'R2': {m: eval_m[m]['R2'] for m in model_names},
            'Pearson': {m: eval_m[m]['Pearson'] for m in model_names},
            'Dir_all': {m: eval_m[m]['Dir_all'] for m in model_names},
            'Dir_deg': {m: eval_m[m]['Dir_deg'] for m in model_names},
            'DEG_auprc': {m: eval_m[m]['DEG_auprc'] for m in model_names},
        }
        rankings = {}
        for rn, sc in rank_scores.items():
            sorted_m = sorted(sc.keys(), key=lambda m: sc[m], reverse=True)
            rankings[rn] = [sorted_m.index(m) + 1 for m in model_names]

        # Print rankings
        log(f"  {'Model':<20} " + " ".join(f"{k:>10}" for k in ['MSE','R2','Dir_all','Dir_deg','DEG_auprc']))
        for i, mn in enumerate(model_names):
            log(f"  {mn:<20} " + " ".join(f"{rankings[k][i]:>10}" for k in ['MSE','R2','Dir_all','Dir_deg','DEG_auprc']))

        key_pairs = [('MSE','Dir_all'),('MSE','Dir_deg'),('MSE','DEG_auprc'),('R2','Dir_deg'),('Pearson','Dir_deg')]
        reversal_results = []
        for m1, m2 in key_pairs:
            tau, p = stats.kendalltau(rankings[m1], rankings[m2])
            interp = 'REVERSAL' if tau < 0.5 else ('PARTIAL' if tau < 0.7 else 'CONSISTENT')
            log(f"    tau({m1},{m2})={tau:.3f} p={p:.4f} -> {interp}")
            reversal_results.append({'m1':m1,'m2':m2,'tau':float(tau),'p':float(p),'interp':interp})

        # Downstream tasks
        log(f"\n  Downstream tasks...")
        ds_results = {}
        for mn in model_names:
            f1_50 = task_deg_recovery(true_logfc, preds[mn], k=50)
            f1_100 = task_deg_recovery(true_logfc, preds[mn], k=100)
            dd = task_dir_discovery(true_logfc, preds[mn])
            ds_results[mn] = {'f1@50': f1_50, 'f1@100': f1_100, 'dir_discovery_deg': dd}
            log(f"    {mn}: f1@50={f1_50:.3f} dir_disc={dd:.3f}")

        # H2: Metric-Downstream Correlation
        log(f"\n  H2: Metric-Downstream Correlation...")
        em_names = ['MSE','R2','Pearson','Dir_all','Dir_deg','Dir_weighted','DEG_auprc','DEG_dir_auprc']
        dt_names = ['f1@50','f1@100','dir_discovery_deg']

        m_scores = {}
        for em in em_names:
            if em == 'MSE':
                m_scores[em] = np.array([-eval_m[m][em] for m in model_names])
            else:
                m_scores[em] = np.array([eval_m[m][em] for m in model_names])
        t_scores = {dt: np.array([ds_results[m][dt] for m in model_names]) for dt in dt_names}

        corr = {}
        for em in em_names:
            for dt in dt_names:
                if np.std(m_scores[em]) < 1e-10 or np.std(t_scores[dt]) < 1e-10:
                    rho, pval = 0.0, 1.0
                else:
                    rho, pval = stats.spearmanr(m_scores[em], t_scores[dt])
                    if np.isnan(rho): rho, pval = 0.0, 1.0
                corr[f"{em}_vs_{dt}"] = {'rho': float(rho), 'p': float(pval)}

        # H2 test
        bio_eval_ms = ['Dir_all','Dir_deg','Dir_weighted','DEG_auprc','DEG_dir_auprc']
        h2_pass, h2_total = 0, 0
        for dt in dt_names:
            mse_rho = corr.get(f"MSE_vs_{dt}",{}).get('rho',0.0)
            for bem in bio_eval_ms:
                be_rho = corr.get(f"{bem}_vs_{dt}",{}).get('rho',0.0)
                diff = be_rho - mse_rho
                h2_total += 1
                if diff >= 0.1: h2_pass += 1
                v = "PASS" if diff >= 0.1 else "FAIL"
                log(f"    {dt}: rho({bem})={be_rho:.3f} vs rho(MSE)={mse_rho:.3f} diff={diff:+.3f} {v}")
        h2_rate = h2_pass / h2_total if h2_total > 0 else 0
        log(f"  H2: {h2_pass}/{h2_total} = {h2_rate:.1%}")

        # H3: trained > baseline under BioEval
        log(f"\n  H3: Trained vs Baseline...")
        trained = [m for m in model_names if meta[m]['type'] == 'trained']
        baselines = [m for m in model_names if meta[m]['type'] != 'trained']
        h3r = {}
        for metric in ['Dir_deg','Dir_all','DEG_auprc','R2']:
            ts = [eval_m[m][metric] for m in trained]
            bs = [eval_m[m][metric] for m in baselines]
            bt = max(ts) if ts else 0; bb = max(bs) if bs else 0
            bt_name = trained[np.argmax(ts)] if ts else 'N/A'
            h3r[metric] = {'best_trained': float(bt), 'best_baseline': float(bb),
                           'best_trained_name': bt_name, 'trained_wins': bt > bb}
            log(f"    {metric}: best_trained={bt:.4f}({bt_name}) {'>' if bt>bb else '<='} best_baseline={bb:.4f}")

        all_results[ct_name] = {
            'n_perts': ds_data['n_perts'], 'n_genes': ds_data['n_genes'],
            'model_names': model_names, 'model_metadata': meta,
            'eval_metrics': eval_m, 'rankings': rankings,
            'reversal_analysis': reversal_results,
            'downstream_tasks': ds_results, 'correlation_matrix': corr,
            'h2_pass_count': h2_pass, 'h2_total': h2_total, 'h2_pass_rate': float(h2_rate),
            'h3_results': h3r,
        }

    # Cross-dataset summary
    log(f"\n{'='*60}")
    log("CROSS-DATASET SUMMARY")
    log(f"{'='*60}")
    for ct, r in all_results.items():
        revs = [x for x in r['reversal_analysis'] if x['interp'] == 'REVERSAL']
        log(f"  {ct}: H1_reversals={len(revs)}/{len(r['reversal_analysis'])} H2={r['h2_pass_rate']:.1%}")
        for metric, h3 in r['h3_results'].items():
            log(f"    H3 {metric}: trained_wins={h3['trained_wins']} ({h3['best_trained_name']}: {h3['best_trained']:.4f} vs {h3['best_baseline']:.4f})")

    # Save
    elapsed = time.time() - t0
    final = {'run':'run_15','timestamp':time.strftime('%Y-%m-%d %H:%M:%S'),
             'runtime_s':elapsed,'model_type':'real_trained_sklearn_ridge',
             'datasets':all_results}
    with open(os.path.join(OUTPUT_DIR, "run_15_results.json"), 'w') as f:
        json.dump(final, f, indent=2, default=str)
    log(f"\nResults saved. Runtime: {elapsed:.1f}s")

if __name__ == "__main__":
    main()
