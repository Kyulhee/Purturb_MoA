"""
run_20d_analyze_results.py
Analyze GEARS results from run_20c (K562+RPE1) and run_20b (Norman)
Produce unified report with H1/H2/H3 assessment including GEARS DL model.

Uses system Python 3.14 (no GEARS import needed — reads .npy files).
Run with: python run_20d_analyze_results.py
"""

import numpy as np
import os
import json
from itertools import combinations

OUT_DIR = "outputs/analysis/run_20"
DATA_DIR = "outputs/analysis/run_04/data/gears_data"

# Import shared metric functions from run_20c
import sys
sys.path.insert(0, OUT_DIR)


def to_array(X):
    return X.toarray() if hasattr(X, 'toarray') else np.asarray(X)


def find_ctrl_key(adata, pert_col):
    vals = list(adata.obs[pert_col].unique())
    for key in ['ctrl', 'control', 'NaN', 'nan', 'ctrl_ess']:
        if key in vals:
            return key
    return adata.obs[pert_col].value_counts().index[0]


def compute_deg_auprc(pred, true, deg_mask):
    if deg_mask.sum() == 0:
        return np.nan
    scores = np.abs(pred)
    labels = deg_mask.astype(float)
    order = np.argsort(-scores)
    ls = labels[order]
    tp = np.cumsum(ls)
    fp = np.cumsum(1 - ls)
    prec = tp / (tp + fp)
    rec = tp / labels.sum()
    prec = np.concatenate([[1], prec])
    rec = np.concatenate([[0], rec])
    return max(0, np.trapezoid(prec, rec))


def compute_f1_at_k(pred, true, deg_mask, k=50):
    if deg_mask.sum() == 0:
        return np.nan
    top_k = np.argsort(-np.abs(pred))[:k]
    pp = np.zeros(len(pred), dtype=bool)
    pp[top_k] = True
    tp = (pp & deg_mask).sum()
    p = tp / k if k > 0 else 0
    r = tp / deg_mask.sum() if deg_mask.sum() > 0 else 0
    return 2 * p * r / (p + r) if (p + r) > 0 else 0


def compute_dir_metrics(pred, true, deg_mask=None):
    ps, ts = np.sign(pred), np.sign(true)
    match = (ps == ts).astype(float)
    match[(ps == 0) & (ts == 0)] = 1.0
    match[(ts == 0) & (ps != 0)] = 0.5
    dir_all = match.mean()
    dir_deg = match[deg_mask].mean() if deg_mask is not None and deg_mask.sum() > 0 else np.nan
    return dir_all, dir_deg


def compute_all_metrics(pred, true, deg_thresh=0.25):
    from scipy import stats
    deg_mask = np.abs(true) > deg_thresh
    mse = np.mean((pred - true) ** 2)
    ss_res = np.sum((pred - true) ** 2)
    ss_tot = np.sum((true - true.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    valid = (np.abs(true) > 1e-10) | (np.abs(pred) > 1e-10)
    pearson = stats.pearsonr(pred[valid], true[valid])[0] if valid.sum() > 2 else 0
    if np.isnan(pearson):
        pearson = 0
    da, dd = compute_dir_metrics(pred, true, deg_mask)
    return {
        'MSE': mse, 'R2': r2, 'Pearson': pearson,
        'Dir_all': da, 'Dir_deg': dd,
        'DEG_auprc': compute_deg_auprc(pred, true, deg_mask),
        'f1@50': compute_f1_at_k(pred, true, deg_mask, k=50),
    }


def ridge_loo_norman(Y_true, perts, ctrl_key='ctrl'):
    """Norman Ridge LOO with binary perturbation features"""
    single_kos = sorted(set(g for p in perts for g in p.split('+') if g != ctrl_key))
    X = np.zeros((len(perts), len(single_kos)))
    for i, p in enumerate(perts):
        for g in p.split('+'):
            if g in single_kos:
                X[i, single_kos.index(g)] = 1.0

    alpha = 1.0
    n, p_feat = X.shape
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    d = s / (s ** 2 + alpha)
    H = (U * d) @ U.T
    beta_hat = Vt.T @ np.diag(d) @ U.T @ Y_true
    Y_hat = X @ beta_hat
    diag_H = np.diag(H)
    mask = (1 - diag_H) > 1e-10
    Y_pred_ridge = np.copy(Y_hat)
    h_corr = (diag_H / (1 - diag_H))[:, None]
    Y_pred_ridge[mask] = Y_hat[mask] + (Y_true[mask] - Y_hat[mask]) * h_corr[mask]
    return Y_pred_ridge


def ridge_loo_replogle(adata, perts, ctrl_key):
    """Replogle Ridge LOO with PCA features"""
    from sklearn.decomposition import PCA

    ctrl_data = to_array(adata[adata.obs['condition'] == ctrl_key].X)
    n_pcs = min(30, ctrl_data.shape[1], ctrl_data.shape[0])
    pca = PCA(n_components=n_pcs).fit(ctrl_data)
    X = np.zeros((len(perts), 2 * n_pcs + 1))
    for i, p in enumerate(perts):
        pdata = to_array(adata[adata.obs['condition'] == p].X)
        pc = pca.transform(pdata)
        X[i] = np.concatenate([pc.mean(0), pc.var(0), [np.log1p(pdata.shape[0])]])

    Y_true = np.array([
        to_array(adata[adata.obs['condition'] == p].X).mean(axis=0)
        - to_array(adata[adata.obs['condition'] == ctrl_key].X).mean(axis=0)
        for p in perts
    ])

    alpha = 1.0
    n, p_feat = X.shape
    if n > p_feat:
        XtX = X.T @ X + alpha * np.eye(p_feat)
        XtX_inv = np.linalg.inv(XtX)
        H = X @ XtX_inv @ X.T
        beta_hat = XtX_inv @ X.T @ Y_true
    else:
        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        d = s / (s ** 2 + alpha)
        H = (U * d) @ U.T
        beta_hat = Vt.T @ np.diag(d) @ U.T @ Y_true

    Y_hat = X @ beta_hat
    diag_H = np.diag(H)
    mask = (1 - diag_H) > 1e-10
    Y_pred_ridge = np.copy(Y_hat)
    h_corr = (diag_H / (1 - diag_H))[:, None]
    Y_pred_ridge[mask] = Y_hat[mask] + (Y_true[mask] - Y_hat[mask]) * h_corr[mask]
    return Y_pred_ridge, Y_true


def compute_kendall_tau(rank1, rank2):
    """Compute Kendall's tau between two rank dictionaries"""
    keys = sorted(set(rank1.keys()) & set(rank2.keys()))
    if len(keys) < 2:
        return np.nan
    r1 = [rank1[k] for k in keys]
    r2 = [rank2[k] for k in keys]
    from scipy.stats import kendalltau
    tau, pval = kendalltau(r1, r2)
    return tau, pval


def main():
    import scanpy as sc

    report_lines = []
    report_lines.append("# BioEval Run 20d — GEARS DL Model Analysis\n")
    report_lines.append(f"**Date**: 2026-05-03\n")

    all_dataset_results = {}

    datasets = [
        ("K562", "replogle_k562_essential"),
        ("RPE1", "replogle_rpe1_essential"),
        ("Norman", "norman"),
    ]

    for cell_type, subdir in datasets:
        print(f"\n{'='*60}")
        print(f"Analyzing {cell_type}...")

        data_path = os.path.join(DATA_DIR, subdir)
        h5ad_path = os.path.join(data_path, "perturb_processed.h5ad")

        if not os.path.exists(h5ad_path):
            print(f"  [SKIP] {h5ad_path} not found")
            continue

        # Load ground truth
        adata = sc.read_h5ad(h5ad_path)
        pert_col = 'condition'
        ctrl_key = find_ctrl_key(adata, pert_col)
        perts = sorted([p for p in adata.obs[pert_col].unique() if p != ctrl_key])
        ctrl_mean = to_array(adata[adata.obs[pert_col] == ctrl_key].X).mean(axis=0)

        if cell_type == "Norman":
            Y_true = np.array([
                to_array(adata[adata.obs[pert_col] == p].X).mean(axis=0) - ctrl_mean
                for p in perts
            ])
        else:
            Y_true = np.array([
                to_array(adata[adata.obs[pert_col] == p].X).mean(axis=0) - ctrl_mean
                for p in perts
            ])

        print(f"  {len(perts)} perturbations, {Y_true.shape[1]} genes")

        model_results = {}

        # === GEARS ===
        gears_pred_path = os.path.join(OUT_DIR, f"{cell_type}_gears_predictions.npy")
        gears_true_path = os.path.join(OUT_DIR, f"{cell_type}_true_effects.npy")

        if os.path.exists(gears_pred_path) and os.path.exists(gears_true_path):
            Y_pred_gears = np.load(gears_pred_path)
            Y_true_gears = np.load(gears_true_path)

            # Check if predictions are non-trivial
            non_zero = np.count_nonzero(Y_pred_gears)
            total = Y_pred_gears.size
            print(f"  GEARS predictions: {non_zero}/{total} non-zero ({100*non_zero/total:.1f}%)")

            if non_zero > 0:
                metrics_list = [compute_all_metrics(Y_pred_gears[i], Y_true_gears[i]) for i in range(len(perts))]
                gears_metrics = {k: np.nanmean([m[k] for m in metrics_list]) for k in metrics_list[0]}
                model_results['GEARS'] = gears_metrics
                print(f"  GEARS: R2={gears_metrics['R2']:.3f}, Dir_deg={gears_metrics['Dir_deg']:.3f}, "
                      f"DEG_auprc={gears_metrics['DEG_auprc']:.3f}")
            else:
                print(f"  GEARS: All-zero predictions — SKIPPED")
                model_results['GEARS'] = None
        else:
            print(f"  GEARS: No prediction files found")
            model_results['GEARS'] = None

        # === Ridge ===
        if cell_type == "Norman":
            Y_pred_ridge = ridge_loo_norman(Y_true, perts, ctrl_key)
        else:
            Y_pred_ridge, Y_true_r = ridge_loo_replogle(adata, perts, ctrl_key)
            Y_true = Y_true_r  # Use the one computed with consistent ctrl_mean

        ridge_metrics_list = [compute_all_metrics(Y_pred_ridge[i], Y_true[i]) for i in range(len(perts))]
        ridge_metrics = {k: np.nanmean([m[k] for m in ridge_metrics_list]) for k in ridge_metrics_list[0]}
        model_results['Ridge'] = ridge_metrics

        # === mean_predictor ===
        mean_pred_ml = [compute_all_metrics(np.zeros(Y_true.shape[1]), Y_true[i]) for i in range(len(perts))]
        mean_pred_m = {k: np.nanmean([m[k] for m in mean_pred_ml]) for k in mean_pred_ml[0]}
        model_results['mean_predictor'] = mean_pred_m

        # === mean_effect ===
        Y_mean_effect = np.tile(Y_true.mean(0, keepdims=True), (Y_true.shape[0], 1))
        mean_eff_ml = [compute_all_metrics(Y_mean_effect[i], Y_true[i]) for i in range(len(perts))]
        mean_eff_m = {k: np.nanmean([m[k] for m in mean_eff_ml]) for k in mean_eff_ml[0]}
        model_results['mean_effect'] = mean_eff_m

        # Print table
        report_lines.append(f"## {cell_type}\n")
        report_lines.append(f"**Perturbations**: {len(perts)}, **Genes**: {Y_true.shape[1]}\n")
        report_lines.append("| Model | R2 | Pearson | Dir_deg | DEG_auprc | f1@50 | MSE |")
        report_lines.append("|-------|----|---------|---------|-----------|-------|-----|")
        for model_name, metrics in model_results.items():
            if metrics is None:
                report_lines.append(f"| {model_name} | N/A | | | | | |")
            else:
                report_lines.append(f"| {model_name} | {metrics['R2']:.3f} | {metrics['Pearson']:.3f} | "
                                   f"{metrics['Dir_deg']:.3f} | {metrics['DEG_auprc']:.3f} | "
                                   f"{metrics['f1@50']:.3f} | {metrics['MSE']:.4f} |")
        report_lines.append("")

        # H3: GEARS vs baselines
        if model_results.get('GEARS') is not None:
            report_lines.append(f"### H3: GEARS vs Baselines ({cell_type})\n")
            for bl_name in ['Ridge', 'mean_predictor', 'mean_effect']:
                bl = model_results[bl_name]
                gears = model_results['GEARS']
                wins = sum(1 for k in ['R2', 'Dir_deg', 'DEG_auprc', 'f1@50'] if gears.get(k, 0) > bl.get(k, 0))
                losses = 4 - wins
                report_lines.append(f"- vs {bl_name}: GEARS wins {wins}/4, loses {losses}/4")
            report_lines.append("")

        # H1: Rank reversal (MSE vs Dir_deg)
        report_lines.append(f"### H1: MSE vs Dir_deg Ranking ({cell_type})\n")
        valid_models = {k: v for k, v in model_results.items() if v is not None}
        if len(valid_models) >= 3:
            mse_rank = {k: i+1 for i, (k, _) in enumerate(
                sorted(valid_models.items(), key=lambda x: x[1]['MSE']))}
            dir_rank = {k: i+1 for i, (k, _) in enumerate(
                sorted(valid_models.items(), key=lambda x: -x[1]['Dir_deg']))}

            report_lines.append(f"- MSE rank: {dict(sorted(mse_rank.items()))}")
            report_lines.append(f"- Dir_deg rank: {dict(sorted(dir_rank.items()))}")

            tau, pval = compute_kendall_tau(mse_rank, dir_rank)
            report_lines.append(f"- Kendall τ(MSE_rank, Dir_deg_rank) = {tau:.3f}, p = {pval:.3f}")
            if abs(tau) < 0.7:
                report_lines.append(f"- **RANK REVERSAL DETECTED** (|τ| < 0.7)")
            else:
                report_lines.append(f"- No significant rank reversal (|τ| >= 0.7)")
        report_lines.append("")

        all_dataset_results[cell_type] = model_results

    # === Summary ===
    report_lines.append("---\n")
    report_lines.append("## Overall H3 Summary\n")

    gears_available = {ct: res for ct, res in all_dataset_results.items()
                       if res.get('GEARS') is not None}

    if gears_available:
        total_metrics = 0
        gears_wins = 0
        for ct, res in gears_available.items():
            for bl_name in ['Ridge', 'mean_predictor', 'mean_effect']:
                if bl_name in res:
                    for k in ['R2', 'Dir_deg', 'DEG_auprc', 'f1@50']:
                        total_metrics += 1
                        if res['GEARS'].get(k, 0) > res[bl_name].get(k, 0):
                            gears_wins += 1

        report_lines.append(f"- GEARS wins {gears_wins}/{total_metrics} comparisons across all datasets and baselines")
        report_lines.append(f"- GEARS datasets available: {', '.join(gears_available.keys())}")
    else:
        report_lines.append("- GEARS results not available for any dataset")

    # Ridge vs mean_predictor (already confirmed from run_16-19)
    ridge_vs_mean = 0
    ridge_total = 0
    for ct, res in all_dataset_results.items():
        if 'Ridge' in res and 'mean_predictor' in res:
            for k in ['R2', 'Dir_deg', 'DEG_auprc', 'f1@50']:
                ridge_total += 1
                if res['Ridge'].get(k, 0) > res['mean_predictor'].get(k, 0):
                    ridge_vs_mean += 1

    report_lines.append(f"- Ridge wins {ridge_vs_mean}/{ridge_total} vs mean_predictor (confirmed from run_16-19)")

    report = '\n'.join(report_lines)
    report_path = os.path.join(OUT_DIR, "run_20d_analysis.md")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)

    print(f"\nReport saved to {report_path}")
    print(f"\n{'='*60}")
    print("SUMMARY:")
    print(report)


if __name__ == '__main__':
    main()
