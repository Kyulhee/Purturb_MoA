"""
BioEval Phase 4 — Downstream Task Correlation Analysis (Run 14)
================================================================
Tests H2: BioEval metrics better predict downstream biological utility than MSE/R2

Downstream tasks:
  T1: DEG Recovery — Precision@k, Recall@k for top-k predicted DEGs
  T2: Hit Prioritization — Spearman rho between predicted and true gene rankings per perturbation
  T3: Effect-Size Recovery — Pearson rho between predicted and true |logFC| for DEGs
  T4: Direction-Guided Discovery — Fraction of true DEGs recovered when following predicted direction

Then:
  Spearman rho(metric, downstream_task) for each metric × each task
  H2: rho(BioEval) > rho(MSE) by >= 0.1
"""

import json
import os
import warnings
import time
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from sklearn.metrics import mean_squared_error, r2_score, average_precision_score
from collections import defaultdict

warnings.filterwarnings('ignore')

OUTPUT_DIR = r"C:\test\llm_project\nexus-science-win\outputs\analysis\run_14"
DATA_DIR = r"C:\test\llm_project\nexus-science-win\outputs\analysis\run_04\data"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# Downstream Task Implementations
# ============================================================

def task_deg_recovery(true_mat, pred_mat, deg_threshold=0.25, k_values=[10, 50, 100, 200]):
    """
    T1: DEG Recovery
    For each perturbation, rank genes by |predicted logFC|, get top-k,
    measure precision and recall against true DEGs.
    """
    n_perts, n_genes = true_mat.shape
    true_deg = np.abs(true_mat) > deg_threshold

    results = {}
    for k in k_values:
        if k > n_genes:
            continue
        precisions = []
        recalls = []
        f1s = []
        for p in range(n_perts):
            # Rank genes by predicted effect size
            pred_rank = np.argsort(np.abs(pred_mat[p]))[::-1]
            top_k_pred = pred_rank[:k]

            # True DEGs for this perturbation
            true_deg_set = set(np.where(true_deg[p])[0])
            if len(true_deg_set) == 0:
                continue

            pred_deg_set = set(top_k_pred)
            tp = len(pred_deg_set & true_deg_set)
            precision = tp / k
            recall = tp / len(true_deg_set) if len(true_deg_set) > 0 else 0.0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

            precisions.append(precision)
            recalls.append(recall)
            f1s.append(f1)

        results[f'precision@{k}'] = float(np.mean(precisions)) if precisions else 0.0
        results[f'recall@{k}'] = float(np.mean(recalls)) if recalls else 0.0
        results[f'f1@{k}'] = float(np.mean(f1s)) if f1s else 0.0

    return results


def task_hit_prioritization(true_mat, pred_mat, deg_threshold=0.25):
    """
    T2: Hit Prioritization
    For each perturbation, compute Spearman rho between predicted and true
    gene rankings (by absolute effect size).
    """
    n_perts, n_genes = true_mat.shape
    rho_values = []
    p_values = []

    for p in range(n_perts):
        true_abs = np.abs(true_mat[p])
        pred_abs = np.abs(pred_mat[p])
        # Need variation in both
        if np.std(true_abs) < 1e-10 or np.std(pred_abs) < 1e-10:
            continue
        rho, pval = stats.spearmanr(true_abs, pred_abs)
        if not np.isnan(rho):
            rho_values.append(rho)
            p_values.append(pval)

    # Also compute for DEGs only
    true_deg = np.abs(true_mat) > deg_threshold
    rho_deg_values = []
    for p in range(n_perts):
        deg_idx = np.where(true_deg[p])[0]
        if len(deg_idx) < 5:
            continue
        true_abs_deg = np.abs(true_mat[p, deg_idx])
        pred_abs_deg = np.abs(pred_mat[p, deg_idx])
        if np.std(true_abs_deg) < 1e-10 or np.std(pred_abs_deg) < 1e-10:
            continue
        rho, _ = stats.spearmanr(true_abs_deg, pred_abs_deg)
        if not np.isnan(rho):
            rho_deg_values.append(rho)

    return {
        'spearman_rho_all_mean': float(np.mean(rho_values)) if rho_values else 0.0,
        'spearman_rho_all_std': float(np.std(rho_values)) if rho_values else 0.0,
        'spearman_rho_deg_mean': float(np.mean(rho_deg_values)) if rho_deg_values else 0.0,
        'spearman_rho_deg_std': float(np.std(rho_deg_values)) if rho_deg_values else 0.0,
        'n_perts_evaluated': len(rho_values),
        'n_perts_deg_evaluated': len(rho_deg_values),
    }


def task_effect_size_recovery(true_mat, pred_mat, deg_threshold=0.25):
    """
    T3: Effect-Size Recovery
    Pearson correlation between predicted and true |logFC| for DEGs.
    This tests whether the model gets the magnitude right on the genes that matter.
    """
    true_deg = np.abs(true_mat) > deg_threshold

    # All genes
    rho_all, _ = stats.pearsonr(np.abs(true_mat.flatten()), np.abs(pred_mat.flatten()))

    # DEGs only
    deg_true = np.abs(true_mat[true_deg])
    deg_pred = np.abs(pred_mat[true_deg])
    if len(deg_true) > 2 and np.std(deg_true) > 1e-10 and np.std(deg_pred) > 1e-10:
        rho_deg, _ = stats.pearsonr(deg_true, deg_pred)
    else:
        rho_deg = np.nan

    # Mean absolute error on DEGs
    mae_deg = float(np.mean(np.abs(true_mat[true_deg] - pred_mat[true_deg]))) if true_deg.sum() > 0 else np.nan

    return {
        'pearson_abs_all': float(rho_all),
        'pearson_abs_deg': float(rho_deg) if not np.isnan(rho_deg) else None,
        'mae_deg': float(mae_deg) if not np.isnan(mae_deg) else None,
    }


def task_direction_guided_discovery(true_mat, pred_mat, deg_threshold=0.25):
    """
    T4: Direction-Guided Discovery
    If you follow the model's predicted direction for each gene,
    what fraction of true DEG directions do you recover?
    Also: for genes where the model predicts UP, how often is it actually UP among true DEGs?
    """
    true_deg = np.abs(true_mat) > deg_threshold

    # Overall directional accuracy on DEGs (same as Dir_deg but per-perturbation variance)
    sign_match = (np.sign(pred_mat) == np.sign(true_mat)).astype(float)
    # Handle zeros
    zero_true = (true_mat == 0)
    zero_pred = (pred_mat == 0)
    sign_match[zero_true & zero_pred] = 1.0
    sign_match[zero_true & ~zero_pred] = 0.5

    # Per-perturbation direction accuracy on DEGs
    dir_acc_per_pert = []
    for p in range(true_mat.shape[0]):
        deg_idx = np.where(true_deg[p])[0]
        if len(deg_idx) == 0:
            continue
        dir_acc_per_pert.append(sign_match[p, deg_idx].mean())

    # "Discovery" scenario: model predicts top-k genes going UP/DOWN,
    # how many are truly UP/DOWN?
    # Use confidence-weighted: rank by |predicted| * sign(predicted)
    discovery_at_50 = []
    discovery_at_100 = []

    for p in range(true_mat.shape[0]):
        deg_idx = np.where(true_deg[p])[0]
        if len(deg_idx) < 3:
            continue

        # Rank by |predicted| among DEGs
        pred_abs_deg = np.abs(pred_mat[p, deg_idx])
        true_signs = np.sign(true_mat[p, deg_idx])
        pred_signs = np.sign(pred_mat[p, deg_idx])

        # Sort by predicted confidence (descending |pred|)
        rank_idx = np.argsort(pred_abs_deg)[::-1]

        for k, store in [(50, discovery_at_50), (100, discovery_at_100)]:
            actual_k = min(k, len(deg_idx))
            if actual_k < 3:
                continue
            top_idx = rank_idx[:actual_k]
            # Fraction where predicted direction matches true direction
            dir_match = (true_signs[top_idx] == pred_signs[top_idx]).mean()
            store.append(dir_match)

    return {
        'dir_acc_deg_mean': float(np.mean(dir_acc_per_pert)) if dir_acc_per_pert else 0.0,
        'dir_acc_deg_std': float(np.std(dir_acc_per_pert)) if dir_acc_per_pert else 0.0,
        'discovery_dir@50': float(np.mean(discovery_at_50)) if discovery_at_50 else 0.0,
        'discovery_dir@100': float(np.mean(discovery_at_100)) if discovery_at_100 else 0.0,
        'n_perts_evaluated': len(dir_acc_per_pert),
    }


# ============================================================
# BioEval metrics (from run_13)
# ============================================================

def bioeval_dir_from_matrices(true_mat, pred_mat, deg_threshold=0.25):
    n_perts, n_genes = true_mat.shape
    sign_true = np.sign(true_mat)
    sign_pred = np.sign(pred_mat)
    sign_match = (sign_true == sign_pred).astype(float)
    zero_true = (true_mat == 0)
    zero_pred = (pred_mat == 0)
    sign_match[zero_true & zero_pred] = 1.0
    sign_match[zero_true & ~zero_pred] = 0.5
    dir_accuracy_all = sign_match.mean()
    deg_mask = np.abs(true_mat) > deg_threshold
    dir_accuracy_deg = sign_match[deg_mask].mean() if deg_mask.sum() > 0 else np.nan
    weights = np.abs(true_mat)
    weights = weights / (weights.sum() + 1e-8)
    dir_accuracy_weighted = (sign_match * weights).sum()
    return {
        'dir_accuracy_all': float(dir_accuracy_all),
        'dir_accuracy_deg': float(dir_accuracy_deg) if not np.isnan(dir_accuracy_deg) else None,
        'dir_accuracy_weighted': float(dir_accuracy_weighted),
    }


def bioeval_cal(true_mat, pred_mat):
    n_perts, n_genes = true_mat.shape
    # Global slope
    t_flat = true_mat.flatten()
    p_flat = pred_mat.flatten()
    if np.std(t_flat) > 1e-10 and np.std(p_flat) > 1e-10:
        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(t_flat, p_flat)
        except ValueError:
            slope, intercept = 0.0, 0.0
    else:
        slope, intercept = 0.0, 0.0

    # Per-gene calibration
    cal_per_gene = []
    for g in range(n_genes):
        t_col = true_mat[:, g]
        p_col = pred_mat[:, g]
        if np.std(t_col) > 1e-10 and np.std(p_col) > 1e-10 and len(np.unique(t_col)) > 1:
            try:
                sg, _, _, _, _ = stats.linregress(t_col, p_col)
                cal_per_gene.append(sg)
            except ValueError:
                cal_per_gene.append(np.nan)
        else:
            cal_per_gene.append(np.nan)

    cal_per_gene = np.array(cal_per_gene)
    valid = ~np.isnan(cal_per_gene)
    cal_slope_dev = float(np.abs(cal_per_gene[valid] - 1.0).mean()) if valid.sum() > 0 else np.nan
    return {
        'global_slope': float(slope),
        'cal_slope_dev': float(cal_slope_dev) if not np.isnan(cal_slope_dev) else None,
    }


def bioeval_deg(true_mat, pred_mat, deg_threshold=0.25):
    y_true_bin = (np.abs(true_mat) > deg_threshold).astype(int).flatten()
    y_pred_score = np.abs(pred_mat).flatten()
    if y_true_bin.sum() > 0 and y_true_bin.sum() < len(y_true_bin):
        auprc = average_precision_score(y_true_bin, y_pred_score)
    else:
        auprc = np.nan
    sign_match = (np.sign(pred_mat) == np.sign(true_mat)).astype(float)
    dir_aware_score = (np.abs(pred_mat) * sign_match).flatten()
    if y_true_bin.sum() > 0 and y_true_bin.sum() < len(y_true_bin):
        dir_aware_auprc = average_precision_score(y_true_bin, dir_aware_score)
    else:
        dir_aware_auprc = np.nan
    return {
        'auprc': float(auprc) if not np.isnan(auprc) else None,
        'dir_aware_auprc': float(dir_aware_auprc) if not np.isnan(dir_aware_auprc) else None,
    }


def compute_existing_metrics(true_mat, pred_mat):
    mse_per_pert = []
    r2_per_pert = []
    pearson_per_pert = []
    for p in range(true_mat.shape[0]):
        y_t = true_mat[p]
        y_p = pred_mat[p]
        mse_per_pert.append(mean_squared_error(y_t, y_p))
        r2_per_pert.append(r2_score(y_t, y_p))
        if np.std(y_t) > 0 and np.std(y_p) > 0:
            r, _ = stats.pearsonr(y_t, y_p)
            pearson_per_pert.append(r)
        else:
            pearson_per_pert.append(np.nan)
    return {
        'mse_mean': float(np.mean(mse_per_pert)),
        'r2_mean': float(np.mean(r2_per_pert)),
        'pearson_mean': float(np.nanmean(pearson_per_pert)),
    }


# ============================================================
# Data Loading + Model Predictions (same as run_13)
# ============================================================

def load_and_prepare(adata, cell_type_name):
    pert_col = None
    for col in ['condition', 'perturbation', 'perturb', 'pert_name']:
        if col in adata.obs.columns:
            pert_col = col
            break
    if pert_col is None:
        print(f"  [{cell_type_name}] No perturbation column found!")
        return None

    pert_values = list(adata.obs[pert_col].unique())
    ctrl_key = None
    for key in ['ctrl', 'control', 'NaN', 'nan', 'ctrl_ess']:
        if key in pert_values:
            ctrl_key = key
            break
    if ctrl_key is None:
        ctrl_key = adata.obs[pert_col].value_counts().index[0]

    print(f"  [{cell_type_name}] Control: '{ctrl_key}', pert_col: '{pert_col}', {len(pert_values)} unique values")

    ctrl_mask = adata.obs[pert_col] == ctrl_key
    ctrl_mean = np.array(adata[ctrl_mask].X.mean(axis=0)).flatten()

    perturbations = [p for p in pert_values if p != ctrl_key]
    actual_means = {}
    for pert in perturbations:
        mask = adata.obs[pert_col] == pert
        if mask.sum() > 0:
            actual_means[pert] = np.array(adata[mask].X.mean(axis=0)).flatten()

    perts_with_data = sorted(actual_means.keys())
    true_logfc = np.array([actual_means[p] - ctrl_mean for p in perts_with_data])

    print(f"  [{cell_type_name}] {len(perts_with_data)} perturbations, {true_logfc.shape[1]} genes")
    print(f"  [{cell_type_name}] logFC: mean={true_logfc.mean():.4f}, std={true_logfc.std():.4f}, DEG_frac={((np.abs(true_logfc) > 0.25).mean()):.4f}")

    return {
        'cell_type': cell_type_name,
        'perturbations': perts_with_data,
        'true_logfc': true_logfc,
        'n_perts': len(perts_with_data),
        'n_genes': true_logfc.shape[1],
    }


def generate_model_predictions(true_logfc):
    gene_std = true_logfc.std(axis=0)
    median_std = np.median(gene_std[gene_std > 0])

    def gene_noise(scale, shape):
        base = np.random.normal(0, 1, shape)
        scale_vec = gene_std.reshape(1, -1) / (median_std + 1e-8)
        return base * scale_vec * scale

    np.random.seed(42)
    model_predictions = {}

    model_predictions['mean_predictor'] = np.zeros_like(true_logfc)
    model_predictions['additive_linear'] = true_logfc * 0.55 + gene_noise(0.03, true_logfc.shape)
    model_predictions['cpa_like'] = true_logfc * 0.65 + gene_noise(0.06, true_logfc.shape)
    model_predictions['gears_like'] = true_logfc * 0.75 + gene_noise(0.08, true_logfc.shape)
    model_predictions['scgpt_like'] = true_logfc * 0.90 + gene_noise(0.04, true_logfc.shape)
    model_predictions['calibrated_noisy'] = true_logfc * 1.0 + gene_noise(0.12, true_logfc.shape)
    model_predictions['over_predictor'] = true_logfc * 1.3 + gene_noise(0.05, true_logfc.shape)

    mean_effect = true_logfc.mean(axis=1, keepdims=True)
    model_predictions['mean_effect_trap'] = mean_effect * 0.3 + gene_noise(0.02, true_logfc.shape)

    deg_mask = np.abs(true_logfc) > 0.25
    partial_flip = true_logfc.copy()
    flip_mask = deg_mask & (np.random.random(true_logfc.shape) < 0.20)
    partial_flip[flip_mask] *= -1
    model_predictions['partial_flip'] = partial_flip * 0.8 + gene_noise(0.04, true_logfc.shape)

    shuffled = true_logfc.copy()
    deg_mask2 = np.abs(true_logfc) > 0.25
    random_sign = np.where(np.random.random(true_logfc.shape) > 0.5, 1.0, -1.0)
    shuffled[deg_mask2] = np.abs(shuffled[deg_mask2]) * random_sign[deg_mask2]
    model_predictions['shuffled_dir'] = shuffled * 0.6 + gene_noise(0.05, true_logfc.shape)

    model_predictions['slight_above_mean'] = true_logfc * 0.15 + gene_noise(0.03, true_logfc.shape)

    return model_predictions


# ============================================================
# Main Phase 4 Pipeline
# ============================================================

def main():
    start_time = time.time()

    print("\n" + "=" * 60)
    print("PHASE 4: Downstream Task Correlation Analysis")
    print("=" * 60)

    # Load data
    print("\nLoading datasets...")
    k562_path = os.path.join(DATA_DIR, "gears_data", "replogle_k562_essential", "perturb_processed.h5ad")
    rpe1_path = os.path.join(DATA_DIR, "gears_data", "replogle_rpe1_essential", "perturb_processed.h5ad")
    norman_path = os.path.join(DATA_DIR, "gears_data", "norman", "perturb_processed.h5ad")

    datasets = {}
    for ct, path in [("K562", k562_path), ("RPE1", rpe1_path), ("Norman", norman_path)]:
        if os.path.exists(path):
            print(f"\nLoading {ct}...")
            adata = sc.read_h5ad(path)
            print(f"  {ct}: {adata.shape[0]} cells, {adata.shape[1]} genes")
            result = load_and_prepare(adata, ct)
            if result:
                datasets[ct] = result
        else:
            print(f"  {ct}: NOT FOUND at {path}")

    if 'K562' not in datasets:
        print("ERROR: K562 data not found.")
        return

    # ============================================================
    # For each dataset: compute metrics + downstream tasks
    # ============================================================
    all_results = {}

    for ct_name, ds_data in datasets.items():
        print(f"\n{'=' * 60}")
        print(f"ANALYZING: {ct_name}")
        print(f"{'=' * 60}")

        true_logfc = ds_data['true_logfc']
        model_predictions = generate_model_predictions(true_logfc)
        model_names = list(model_predictions.keys())
        print(f"Models ({len(model_names)}): {model_names}")

        # ---- Compute evaluation metrics (independent vars) ----
        print(f"\n  Computing evaluation metrics...")
        eval_metrics = {}
        for model_name in model_names:
            pred_mat = model_predictions[model_name]

            # Existing metrics
            existing = compute_existing_metrics(true_logfc, pred_mat)

            # BioEval metrics
            dir_result = bioeval_dir_from_matrices(true_logfc, pred_mat, deg_threshold=0.25)
            cal_result = bioeval_cal(true_logfc, pred_mat)
            deg_result = bioeval_deg(true_logfc, pred_mat, deg_threshold=0.25)

            eval_metrics[model_name] = {
                'MSE': existing['mse_mean'],
                'R2': existing['r2_mean'],
                'Pearson': existing['pearson_mean'],
                'Dir_all': dir_result['dir_accuracy_all'],
                'Dir_deg': dir_result['dir_accuracy_deg'] if dir_result['dir_accuracy_deg'] is not None else 0.0,
                'Dir_weighted': dir_result['dir_accuracy_weighted'],
                'Cal_slope_dev': cal_result['cal_slope_dev'] if cal_result['cal_slope_dev'] is not None else 1.0,
                'DEG_auprc': deg_result['auprc'] if deg_result['auprc'] is not None else 0.0,
                'DEG_dir_auprc': deg_result['dir_aware_auprc'] if deg_result['dir_aware_auprc'] is not None else 0.0,
            }

        # ---- Compute downstream tasks (dependent vars) ----
        print(f"  Computing downstream tasks...")
        downstream_results = {}
        for model_name in model_names:
            pred_mat = model_predictions[model_name]

            t1 = task_deg_recovery(true_logfc, pred_mat)
            t2 = task_hit_prioritization(true_logfc, pred_mat)
            t3 = task_effect_size_recovery(true_logfc, pred_mat)
            t4 = task_direction_guided_discovery(true_logfc, pred_mat)

            downstream_results[model_name] = {
                # T1: DEG Recovery
                'prec@50': t1.get('precision@50', 0.0),
                'recall@50': t1.get('recall@50', 0.0),
                'f1@50': t1.get('f1@50', 0.0),
                'prec@100': t1.get('precision@100', 0.0),
                'recall@100': t1.get('recall@100', 0.0),
                'f1@100': t1.get('f1@100', 0.0),
                # T2: Hit Prioritization
                'spearman_rho_all': t2['spearman_rho_all_mean'],
                'spearman_rho_deg': t2['spearman_rho_deg_mean'],
                # T3: Effect-Size Recovery
                'pearson_abs_all': t3['pearson_abs_all'],
                'pearson_abs_deg': t3['pearson_abs_deg'] if t3['pearson_abs_deg'] is not None else 0.0,
                'mae_deg': t3['mae_deg'] if t3['mae_deg'] is not None else 0.0,
                # T4: Direction-Guided Discovery
                'dir_discovery_deg': t4['dir_acc_deg_mean'],
                'discovery_dir@50': t4['discovery_dir@50'],
                'discovery_dir@100': t4['discovery_dir@100'],
            }

            print(f"    {model_name}: prec@50={t1.get('precision@50', 0.0):.3f}, "
                  f"rho_deg={t2['spearman_rho_deg_mean']:.3f}, "
                  f"dir_disc={t4['dir_acc_deg_mean']:.3f}")

        # ---- H2: Metric-Downstream Correlation ----
        print(f"\n  H2: Metric-Downstream Correlation...")

        eval_metric_names = ['MSE', 'R2', 'Pearson', 'Dir_all', 'Dir_deg', 'Dir_weighted',
                             'Cal_slope_dev', 'DEG_auprc', 'DEG_dir_auprc']

        downstream_task_names = ['prec@50', 'recall@50', 'f1@50',
                                'prec@100', 'recall@100', 'f1@100',
                                'spearman_rho_all', 'spearman_rho_deg',
                                'pearson_abs_all', 'pearson_abs_deg',
                                'dir_discovery_deg', 'discovery_dir@50', 'discovery_dir@100']

        # Build score vectors
        metric_scores = {}
        for em in eval_metric_names:
            # For MSE and Cal_slope_dev, lower is better -> negate for correlation
            if em in ['MSE', 'Cal_slope_dev']:
                metric_scores[em] = np.array([-eval_metrics[m][em] for m in model_names])
            else:
                metric_scores[em] = np.array([eval_metrics[m][em] for m in model_names])

        task_scores = {}
        for dt in downstream_task_names:
            task_scores[dt] = np.array([downstream_results[m][dt] for m in model_names])

        # Compute Spearman rho between each metric and each downstream task
        correlation_matrix = {}
        h2_results = []

        # Key downstream tasks for H2 testing
        key_downstream = ['f1@50', 'f1@100', 'spearman_rho_deg', 'pearson_abs_deg',
                          'dir_discovery_deg', 'discovery_dir@100']

        print(f"\n  Key Metric-Downstream Correlations (Spearman rho):")
        print(f"  {'Metric':<18}", end="")
        for dt in key_downstream:
            print(f"  {dt:>16}", end="")
        print()

        for em in eval_metric_names:
            row = f"  {em:<18}"
            for dt in key_downstream:
                if np.std(metric_scores[em]) < 1e-10 or np.std(task_scores[dt]) < 1e-10:
                    rho = 0.0
                    pval = 1.0
                else:
                    rho, pval = stats.spearmanr(metric_scores[em], task_scores[dt])
                    if np.isnan(rho):
                        rho = 0.0
                        pval = 1.0
                correlation_matrix[f"{em}_vs_{dt}"] = {'rho': float(rho), 'p': float(pval)}
                row += f"  {rho:>16.3f}"
            print(row)

        # H2 test: BioEval metrics vs MSE
        print(f"\n  H2 Test: rho(BioEval) vs rho(MSE) for each downstream task:")
        bioeval_metrics = ['Dir_all', 'Dir_deg', 'Dir_weighted', 'DEG_auprc', 'DEG_dir_auprc']
        h2_pass_count = 0
        h2_total = 0

        for dt in key_downstream:
            mse_rho = correlation_matrix.get(f"MSE_vs_{dt}", {}).get('rho', 0.0)
            print(f"\n    {dt}: rho(MSE) = {mse_rho:.3f}")
            for bem in bioeval_metrics:
                be_rho = correlation_matrix.get(f"{bem}_vs_{dt}", {}).get('rho', 0.0)
                diff = be_rho - mse_rho
                h2_total += 1
                passed = diff >= 0.1
                if passed:
                    h2_pass_count += 1
                verdict = "PASS" if passed else "FAIL"
                print(f"      rho({bem}) = {be_rho:.3f} (diff = {diff:+.3f}) -> {verdict}")

        h2_pass_rate = h2_pass_count / h2_total if h2_total > 0 else 0
        print(f"\n  H2 Summary: {h2_pass_count}/{h2_total} comparisons pass (rate = {h2_pass_rate:.1%})")
        h2_overall = "SUPPORTED" if h2_pass_rate >= 0.5 else "NOT SUPPORTED"
        print(f"  H2 Overall: {h2_overall} (threshold: >= 50% of comparisons with rho(BioEval) >= rho(MSE) + 0.1)")

        # ---- Per-model downstream task summary table ----
        print(f"\n  Per-Model Downstream Task Summary ({ct_name}):")
        print(f"  {'Model':<22} {'prec@50':>8} {'f1@50':>8} {'rho_deg':>8} {'dir_disc':>8} {'disc@100':>8}")
        for model_name in model_names:
            dr = downstream_results[model_name]
            print(f"  {model_name:<22} {dr['prec@50']:>8.3f} {dr['f1@50']:>8.3f} "
                  f"{dr['spearman_rho_deg']:>8.3f} {dr['dir_discovery_deg']:>8.3f} "
                  f"{dr['discovery_dir@100']:>8.3f}")

        # Store results
        all_results[ct_name] = {
            'eval_metrics': {m: eval_metrics[m] for m in model_names},
            'downstream_tasks': {m: downstream_results[m] for m in model_names},
            'correlation_matrix': correlation_matrix,
            'h2_pass_count': h2_pass_count,
            'h2_total': h2_total,
            'h2_pass_rate': float(h2_pass_rate),
            'h2_overall': h2_overall,
            'model_names': model_names,
        }

    # ============================================================
    # Cross-dataset H2 summary
    # ============================================================
    print(f"\n{'=' * 60}")
    print("CROSS-DATASET H2 SUMMARY")
    print(f"{'=' * 60}")

    for ct_name, ct_results in all_results.items():
        print(f"  {ct_name}: H2 pass rate = {ct_results['h2_pass_rate']:.1%} -> {ct_results['h2_overall']}")

    # Aggregate
    total_pass = sum(r['h2_pass_count'] for r in all_results.values())
    total_all = sum(r['h2_total'] for r in all_results.values())
    overall_rate = total_pass / total_all if total_all > 0 else 0
    print(f"\n  Aggregate: {total_pass}/{total_all} = {overall_rate:.1%}")

    # ============================================================
    # Key insight: MSE vs Dir_deg downstream correlation comparison
    # ============================================================
    print(f"\n{'=' * 60}")
    print("KEY COMPARISON: MSE vs Dir_deg downstream correlation")
    print(f"{'=' * 60}")

    key_tasks = ['f1@50', 'f1@100', 'spearman_rho_deg', 'dir_discovery_deg', 'discovery_dir@100']
    for ct_name, ct_results in all_results.items():
        print(f"\n  {ct_name}:")
        for dt in key_tasks:
            mse_rho = ct_results['correlation_matrix'].get(f"MSE_vs_{dt}", {}).get('rho', 0.0)
            dir_rho = ct_results['correlation_matrix'].get(f"Dir_deg_vs_{dt}", {}).get('rho', 0.0)
            diff = dir_rho - mse_rho
            arrow = ">>>" if diff >= 0.1 else (">" if diff > 0 else "<=")
            print(f"    {dt}: rho(MSE)={mse_rho:.3f}, rho(Dir_deg)={dir_rho:.3f} (diff={diff:+.3f}) {arrow}")

    # ============================================================
    # Save results
    # ============================================================
    elapsed = time.time() - start_time

    final_results = {
        'run': 'run_14_phase4',
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'runtime_s': elapsed,
        'h2_test': {
            'threshold': 'rho(BioEval) >= rho(MSE) + 0.1',
            'pass_rate_threshold': 0.5,
            'aggregate_pass_rate': float(overall_rate),
            'aggregate_pass_count': total_pass,
            'aggregate_total': total_all,
        },
        'datasets': {},
    }

    for ct_name, ct_results in all_results.items():
        final_results['datasets'][ct_name] = {
            'h2_pass_rate': ct_results['h2_pass_rate'],
            'h2_overall': ct_results['h2_overall'],
            'correlation_matrix': ct_results['correlation_matrix'],
            'eval_metrics': ct_results['eval_metrics'],
            'downstream_tasks': ct_results['downstream_tasks'],
        }

    results_path = os.path.join(OUTPUT_DIR, "run_14_results.json")
    with open(results_path, 'w') as f:
        json.dump(final_results, f, indent=2, default=str)
    print(f"\nResults saved to: {results_path}")

    print(f"\nTotal runtime: {elapsed:.1f}s")
    print("Phase 4 analysis complete!")


if __name__ == "__main__":
    main()
