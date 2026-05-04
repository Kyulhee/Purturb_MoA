"""
run_20c_k562_rpe1.py
B2: Train GEARS on K562 + RPE1 only (Norman handled separately via run_20b)

Single-process run to avoid GPU contention.
Saves intermediate results after each dataset completes.

Uses ai_env conda environment (Python 3.11, PyTorch 2.5.1, PyG 2.7.0).
Run with: "C:/Users/hgh97/miniconda3/envs/ai_env/python.exe" run_20c_k562_rpe1.py
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

# Only K562 + RPE1 (Norman via run_20b)
DATASETS = [
    ("K562", "replogle_k562_essential"),
    ("RPE1", "replogle_rpe1_essential"),
]


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
    if torch.cuda.is_available():
        print(f"GPU: {torch.cuda.get_device_name(0)}, VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    from gears import PertData, GEARS

    all_results = {}

    for cell_type, subdir in DATASETS:
        print(f"\n{'='*60}")
        print(f"Processing {cell_type}...")
        t_start = time.time()

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
        Y_true = np.array([
            to_array(adata[adata.obs[pert_col] == p].X).mean(axis=0) - ctrl_mean
            for p in perts
        ])
        print(f"  {len(perts)} perturbations, {Y_true.shape[1]} genes")

        # ============================================================
        # GEARS Training
        # ============================================================
        print(f"  Training GEARS model...")

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

            # Extract predictions — use saved_pred from training, then predict missing ones
            # Key fix: pert format is "gene+ctrl", need to extract just the gene names
            # and exclude 'ctrl' from the perturbation gene list for GEARS predict()
            preds = {}

            # First try saved predictions from training (test set predictions)
            # GEARS predict returns absolute expression, so subtract ctrl_mean to get perturbation effect
            for key, val in gears_model.saved_pred.items():
                preds[key] = np.asarray(val).flatten() - ctrl_mean

            print(f"  Saved predictions from training: {len(preds)}")

            # Predict remaining perturbations
            for pert in perts:
                pert_genes = [g for g in pert.split('+') if g != ctrl_key and g.lower() != 'ctrl' and g.lower() != 'nan']
                if not pert_genes:
                    preds[pert] = None
                    continue
                key = '_'.join(pert_genes)
                if key in preds:
                    preds[pert] = preds[key]
                    continue
                try:
                    pred = gears_model.predict([pert_genes])
                    if key in pred:
                        preds[pert] = np.asarray(pred[key]).flatten() - ctrl_mean
                    else:
                        preds[pert] = None
                except Exception as e:
                    print(f"    [WARN] predict {pert} (genes={pert_genes}): {e}")
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

            print(f"  GEARS predictions: {n_valid}/{len(perts)} valid")

            # Save immediately after each dataset
            np.save(os.path.join(OUT_DIR, f"{cell_type}_gears_predictions.npy"), Y_pred_gears)
            np.save(os.path.join(OUT_DIR, f"{cell_type}_true_effects.npy"), Y_true)
            print(f"  Saved predictions for {cell_type}")

            metrics_list = [compute_all_metrics(Y_pred_gears[i], Y_true[i]) for i in range(len(perts))]
            gears_metrics = {k: np.nanmean([m[k] for m in metrics_list]) for k in metrics_list[0]}

            print(f"  GEARS metrics: R2={gears_metrics['R2']:.3f}, Dir_deg={gears_metrics['Dir_deg']:.3f}, "
                  f"DEG_auprc={gears_metrics['DEG_auprc']:.3f}")

            results = {'GEARS': gears_metrics}

        except Exception as e:
            print(f"  [ERROR] GEARS training failed: {e}")
            import traceback
            traceback.print_exc()
            results = {'GEARS': 'FAILED'}

        # ============================================================
        # Ridge baseline
        # ============================================================
        print(f"  Computing Ridge baseline metrics...")
        try:
            from sklearn.decomposition import PCA

            ctrl_data = to_array(adata[adata.obs[pert_col] == ctrl_key].X)
            n_pcs = min(30, ctrl_data.shape[1], ctrl_data.shape[0])
            pca = PCA(n_components=n_pcs).fit(ctrl_data)
            X = np.zeros((len(perts), 2 * n_pcs + 1))
            for i, p in enumerate(perts):
                pdata = to_array(adata[adata.obs[pert_col] == p].X)
                pc = pca.transform(pdata)
                X[i] = np.concatenate([pc.mean(0), pc.var(0), [np.log1p(pdata.shape[0])]])

            # Ridge LOO
            n, p_feat = X.shape
            alpha = 1.0
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

            ridge_metrics_list = [compute_all_metrics(Y_pred_ridge[i], Y_true[i]) for i in range(len(perts))]
            ridge_metrics = {k: np.nanmean([m[k] for m in ridge_metrics_list]) for k in ridge_metrics_list[0]}
            results['Ridge'] = ridge_metrics

            # mean_predictor
            mean_pred_metrics_list = [compute_all_metrics(np.zeros(Y_true.shape[1]), Y_true[i]) for i in range(len(perts))]
            mean_pred_metrics = {k: np.nanmean([m[k] for m in mean_pred_metrics_list]) for k in mean_pred_metrics_list[0]}
            results['mean_predictor'] = mean_pred_metrics

            # mean_effect
            Y_mean_effect = np.tile(Y_true.mean(0, keepdims=True), (Y_true.shape[0], 1))
            mean_eff_metrics_list = [compute_all_metrics(Y_mean_effect[i], Y_true[i]) for i in range(len(perts))]
            mean_eff_metrics = {k: np.nanmean([m[k] for m in mean_eff_metrics_list]) for k in mean_eff_metrics_list[0]}
            results['mean_effect'] = mean_eff_metrics

        except Exception as e:
            print(f"  [ERROR] Ridge computation failed: {e}")
            import traceback
            traceback.print_exc()

        all_results[cell_type] = results
        elapsed_ds = time.time() - t_start
        print(f"  {cell_type} completed in {elapsed_ds:.1f}s")

        # Save intermediate report after each dataset
        R = [f"# BioEval Run 20c Report — GEARS K562+RPE1\n"]
        R.append(f"**Date**: 2026-05-02")
        R.append(f"**Runtime so far**: {time.time()-t0:.1f}s\n")
        for ct, res in all_results.items():
            R.append(f"## {ct}\n")
            R.append("| Model | R2 | Dir_deg | DEG_auprc | f1@50 | MSE |")
            R.append("|-------|----|---------|-----------|-------|-----|")
            for model_name, metrics in res.items():
                if isinstance(metrics, str):
                    R.append(f"| {model_name} | {metrics} | | | | |")
                else:
                    R.append(f"| {model_name} | {metrics.get('R2', 0):.3f} | {metrics.get('Dir_deg', 0):.3f} | "
                             f"{metrics.get('DEG_auprc', 0):.3f} | {metrics.get('f1@50', 0):.3f} | "
                             f"{metrics.get('MSE', 0):.4f} |")
            R.append("")
        report = '\n'.join(R)
        with open(os.path.join(OUT_DIR, "run_20c_report.md"), 'w', encoding='utf-8') as f:
            f.write(report)

    print(f"\nTotal runtime: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
