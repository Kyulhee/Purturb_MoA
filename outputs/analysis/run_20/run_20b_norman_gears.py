"""
run_20b_norman_gears.py
B2: Train GEARS on Norman dataset only (GO graph CSV skip patch applied)

The main run_20 failed on Norman due to corrupted GO graph tar.gz.
This script re-runs Norman only, leveraging the patched utils.py
that skips tar extraction when go_essential_all.csv already exists.

Uses ai_env conda environment (Python 3.11, PyTorch 2.5.1, PyG 2.7.0).
Run with: "C:/Users/hgh97/miniconda3/envs/ai_env/python.exe" run_20b_norman_gears.py
"""

import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')

DATA_DIR = "outputs/analysis/run_04/data/gears_data"
OUT_DIR = "outputs/analysis/run_20"
os.makedirs(OUT_DIR, exist_ok=True)

GEARS_CONFIG = {
    'device': 'cuda',
    'seed': 42,
    'num_epochs': 5,
    'lr': 1e-3,
    'weight_decay': 1e-5,
    'hidden_size': 64,
    'num_heads': 2,
    'num_layers': 2,
}


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


def main():
    t0 = time.time()
    import torch
    print(f"PyTorch {torch.__version__}, CUDA={torch.cuda.is_available()}")

    from gears import PertData, GEARS

    cell_type = "Norman"
    subdir = "norman"

    data_path = os.path.join(DATA_DIR, subdir)
    h5ad_path = os.path.join(data_path, "perturb_processed.h5ad")

    if not os.path.exists(h5ad_path):
        print(f"[ERROR] {h5ad_path} not found")
        return

    # Load ground truth
    adata = sc.read_h5ad(h5ad_path)
    pert_col = 'condition'
    ctrl_key = find_ctrl_key(adata, pert_col)
    perts = sorted([p for p in adata.obs[pert_col].unique() if p != ctrl_key])
    ctrl_mean = to_array(adata[adata.obs[pert_col] == ctrl_key].X).mean(axis=0)
    Y_true = np.array([
        to_array(adata[adata.obs[pert_col] == p].X).mean(axis=0) - ctrl_mean
        for p in perts
    ])
    print(f"{len(perts)} perturbations, {Y_true.shape[1]} genes")

    # GEARS Training
    print("Training GEARS model on Norman...")

    try:
        pert_data = PertData(DATA_DIR)
        pert_data.load(data_name=subdir)
        pert_data.prepare_split(split='simulation', seed=42)
        pert_data.get_dataloader(batch_size=32)

        gears_model = GEARS(pert_data, device='cuda')
        gears_model.model_initialize(
            hidden_size=GEARS_CONFIG['hidden_size'],
            num_go_gnn_layers=GEARS_CONFIG['num_layers'],
            num_gene_gnn_layers=GEARS_CONFIG['num_layers'],
        )

        gears_model.train(
            epochs=GEARS_CONFIG['num_epochs'],
            lr=GEARS_CONFIG['lr'],
            weight_decay=GEARS_CONFIG['weight_decay'],
        )

        # Extract predictions
        preds = {}
        for pert in perts:
            pert_genes = pert.split('+')
            try:
                pred = gears_model.predict([pert_genes])
                key = '_'.join(pert_genes)
                if key in pred:
                    preds[pert] = np.asarray(pred[key]).flatten() - ctrl_mean
                else:
                    preds[pert] = None
            except Exception as e:
                print(f"  [WARN] predict {pert}: {e}")
                preds[pert] = None

        Y_pred_gears = np.zeros_like(Y_true)
        n_valid = 0
        for i, pert in enumerate(perts):
            if preds.get(pert) is not None:
                p = np.asarray(preds[pert]).flatten()
                if len(p) == Y_true.shape[1]:
                    Y_pred_gears[i] = p
                    n_valid += 1
                else:
                    Y_pred_gears[i] = 0
            else:
                Y_pred_gears[i] = 0

        print(f"GEARS predictions: {n_valid}/{len(perts)} valid")

        metrics_list = [compute_all_metrics(Y_pred_gears[i], Y_true[i]) for i in range(len(perts))]
        gears_metrics = {k: np.nanmean([m[k] for m in metrics_list]) for k in metrics_list[0]}

        print(f"GEARS metrics: R2={gears_metrics['R2']:.3f}, Dir_deg={gears_metrics['Dir_deg']:.3f}, "
              f"DEG_auprc={gears_metrics['DEG_auprc']:.3f}")

        np.save(os.path.join(OUT_DIR, "Norman_gears_predictions.npy"), Y_pred_gears)
        np.save(os.path.join(OUT_DIR, "Norman_true_effects.npy"), Y_true)

        # Save report
        R = []
        R.append("# BioEval Run 20b Report — Norman GEARS (GO graph fix)\n")
        R.append(f"**Date**: 2026-05-02")
        R.append(f"**Runtime**: {time.time()-t0:.1f}s\n")
        R.append("## Norman GEARS Results\n")
        R.append("| Model | R2 | Dir_deg | DEG_auprc | f1@50 | MSE |")
        R.append("|-------|----|---------|-----------|-------|-----|")
        R.append(f"| GEARS | {gears_metrics['R2']:.3f} | {gears_metrics['Dir_deg']:.3f} | "
                 f"{gears_metrics['DEG_auprc']:.3f} | {gears_metrics['f1@50']:.3f} | "
                 f"{gears_metrics['MSE']:.4f} |")

        # Ridge baseline for comparison
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

        ridge_ml = [compute_all_metrics(Y_pred_ridge[i], Y_true[i]) for i in range(len(perts))]
        ridge_m = {k: np.nanmean([m[k] for m in ridge_ml]) for k in ridge_ml[0]}
        R.append(f"| Ridge | {ridge_m['R2']:.3f} | {ridge_m['Dir_deg']:.3f} | "
                 f"{ridge_m['DEG_auprc']:.3f} | {ridge_m['f1@50']:.3f} | "
                 f"{ridge_m['MSE']:.4f} |")

        # mean_predictor
        mean_ml = [compute_all_metrics(np.zeros(Y_true.shape[1]), Y_true[i]) for i in range(len(perts))]
        mean_m = {k: np.nanmean([m[k] for m in mean_ml]) for k in mean_ml[0]}
        R.append(f"| mean_predictor | {mean_m['R2']:.3f} | {mean_m['Dir_deg']:.3f} | "
                 f"{mean_m['DEG_auprc']:.3f} | {mean_m['f1@50']:.3f} | "
                 f"{mean_m['MSE']:.4f} |")

        # H3 comparison
        R.append("\n### H3: GEARS vs Baselines\n")
        for bl_name, bl_m in [('Ridge', ridge_m), ('mean_predictor', mean_m)]:
            wins = sum(1 for k in ['R2', 'Dir_deg', 'DEG_auprc', 'f1@50'] if gears_metrics.get(k, 0) > bl_m.get(k, 0))
            R.append(f"- vs {bl_name}: GEARS wins {wins}/4 metrics")

        report = '\n'.join(R)
        with open(os.path.join(OUT_DIR, "run_20b_report.md"), 'w', encoding='utf-8') as f:
            f.write(report)
        print(f"\nReport saved to {OUT_DIR}/run_20b_report.md")

    except Exception as e:
        print(f"[ERROR] Norman GEARS training failed: {e}")
        import traceback
        traceback.print_exc()

    print(f"Total runtime: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
