"""
BioEval Analysis — Run 13 (v2)
========================
Phase 1: Load data + generate realistic model predictions
Phase 2: Implement BioEval metrics + compute existing metrics
Phase 3: Metric-ranking reversal analysis (RQ1 + RQ3)
Phase 4: Sensitivity S1 + cross-cell-type analysis
"""

import json
import os
import warnings
import time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import stats
from sklearn.metrics import mean_squared_error, r2_score, average_precision_score
from collections import defaultdict

warnings.filterwarnings('ignore')

OUTPUT_DIR = r"C:\test\llm_project\nexus-science-win\outputs\analysis\run_13"
DATA_DIR = r"C:\test\llm_project\nexus-science-win\outputs\analysis\run_04\data"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# Phase 2: BioEval Metrics Implementation
# ============================================================

def bioeval_dir_from_matrices(true_mat, pred_mat, deg_threshold=0.25):
    """
    BioEval-Dir from matrices (n_perturbations, n_genes)
    """
    n_perts, n_genes = true_mat.shape

    # Sign match
    sign_true = np.sign(true_mat)
    sign_pred = np.sign(pred_mat)
    sign_match = (sign_true == sign_pred).astype(float)

    # Handle zeros
    zero_true = (true_mat == 0)
    zero_pred = (pred_mat == 0)
    sign_match[zero_true & zero_pred] = 1.0
    sign_match[zero_true & ~zero_pred] = 0.5

    # All genes
    dir_accuracy_all = sign_match.mean()
    dir_per_pert = sign_match.mean(axis=1)
    dir_per_gene = sign_match.mean(axis=0)

    # DEG-only
    deg_mask = np.abs(true_mat) > deg_threshold
    dir_accuracy_deg = sign_match[deg_mask].mean() if deg_mask.sum() > 0 else np.nan

    # Magnitude ratio
    eps = 1e-8
    ratio = np.abs(pred_mat) / (np.abs(true_mat) + eps)
    ratio = np.clip(ratio, 0, 10)
    mag_ratio_all = np.median(ratio)
    mag_ratio_deg = np.median(ratio[deg_mask]) if deg_mask.sum() > 0 else np.nan

    # Weighted direction accuracy (weight by |true|)
    weights = np.abs(true_mat)
    weights = weights / (weights.sum() + eps)
    dir_accuracy_weighted = (sign_match * weights).sum()

    return {
        'dir_accuracy_all': float(dir_accuracy_all),
        'dir_accuracy_deg': float(dir_accuracy_deg) if not np.isnan(dir_accuracy_deg) else None,
        'dir_accuracy_weighted': float(dir_accuracy_weighted),
        'mag_ratio_all': float(mag_ratio_all),
        'mag_ratio_deg': float(mag_ratio_deg) if not np.isnan(mag_ratio_deg) else None,
        'dir_per_pert_mean': float(dir_per_pert.mean()),
        'dir_per_pert_std': float(dir_per_pert.std()),
        'dir_per_gene_mean': float(dir_per_gene.mean()),
        'dir_per_gene_std': float(dir_per_gene.std()),
        'deg_threshold': deg_threshold,
        'n_deg': int(deg_mask.sum()),
        'frac_deg': float(deg_mask.mean()),
        'dir_per_pert': dir_per_pert.tolist(),
        'dir_per_gene': dir_per_gene.tolist(),
    }

def bioeval_cal(true_mat, pred_mat):
    """
    BioEval-Cal: Effect-size calibration analysis
    """
    n_perts, n_genes = true_mat.shape

    # Global calibration
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        true_mat.flatten(), pred_mat.flatten()
    )

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

    underpredict_frac = float((cal_per_gene[valid] < 0.8).mean())
    overpredict_frac = float((cal_per_gene[valid] > 1.2).mean())
    cal_in_range_frac = float(((cal_per_gene[valid] >= 0.8) & (cal_per_gene[valid] <= 1.2)).mean())

    return {
        'global_slope': float(slope),
        'global_intercept': float(intercept),
        'global_r': float(r_value),
        'global_p': float(p_value),
        'underpredict_frac': underpredict_frac,
        'overpredict_frac': overpredict_frac,
        'cal_in_range_frac': cal_in_range_frac,
        'gene_slope_mean': float(np.nanmean(cal_per_gene)),
        'gene_slope_std': float(np.nanstd(cal_per_gene)),
    }

def bioeval_deg(true_mat, pred_mat, deg_threshold=0.25):
    """
    BioEval-DEG: DEG precision-recall with direction awareness
    """
    y_true_bin = (np.abs(true_mat) > deg_threshold).astype(int).flatten()

    # Standard AUPRC
    y_pred_score = np.abs(pred_mat).flatten()
    if y_true_bin.sum() > 0 and y_true_bin.sum() < len(y_true_bin):
        auprc = average_precision_score(y_true_bin, y_pred_score)
    else:
        auprc = np.nan

    # Direction-aware
    sign_match = (np.sign(pred_mat) == np.sign(true_mat)).astype(float)
    dir_aware_score = (np.abs(pred_mat) * sign_match).flatten()
    if y_true_bin.sum() > 0 and y_true_bin.sum() < len(y_true_bin):
        dir_aware_auprc = average_precision_score(y_true_bin, dir_aware_score)
    else:
        dir_aware_auprc = np.nan

    # DEG overlap (Jaccard)
    pred_deg = np.abs(pred_mat) > deg_threshold
    obs_deg = np.abs(true_mat) > deg_threshold
    intersection = (pred_deg & obs_deg).sum()
    union = (pred_deg | obs_deg).sum()
    jaccard = intersection / union if union > 0 else 0.0

    # Direction-aware DEG overlap
    dir_match_deg = pred_deg & obs_deg & (np.sign(pred_mat) == np.sign(true_mat))
    dir_jaccard = dir_match_deg.sum() / union if union > 0 else 0.0

    return {
        'auprc': float(auprc) if not np.isnan(auprc) else None,
        'dir_aware_auprc': float(dir_aware_auprc) if not np.isnan(dir_aware_auprc) else None,
        'jaccard': float(jaccard),
        'dir_jaccard': float(dir_jaccard),
        'deg_threshold': deg_threshold,
    }

def compute_existing_metrics(true_mat, pred_mat):
    """Compute MSE, R2, Pearson"""
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
        'mse_std': float(np.std(mse_per_pert)),
        'r2_mean': float(np.mean(r2_per_pert)),
        'r2_std': float(np.std(r2_per_pert)),
        'pearson_mean': float(np.nanmean(pearson_per_pert)),
        'pearson_std': float(np.nanstd(pearson_per_pert)),
    }

# ============================================================
# Phase 3: Metric-Ranking Reversal Analysis
# ============================================================

def kendall_rank_analysis(rankings_dict):
    """Compute Kendall tau between all pairs of metric rankings."""
    metric_names = list(rankings_dict.keys())
    n_metrics = len(metric_names)
    results = []
    for i in range(n_metrics):
        for j in range(i+1, n_metrics):
            m1, m2 = metric_names[i], metric_names[j]
            r1 = rankings_dict[m1]
            r2 = rankings_dict[m2]
            tau, p = stats.kendalltau(r1, r2)
            results.append({
                'metric_1': m1,
                'metric_2': m2,
                'kendall_tau': float(tau),
                'p_value': float(p),
                'interpretation': (
                    'REVERSAL' if tau < 0.5 else
                    ('PARTIAL' if tau < 0.7 else 'CONSISTENT')
                )
            })
    return results

# ============================================================
# Data loading helper
# ============================================================

def load_and_prepare(adata, cell_type_name):
    """Load adata, find control, compute logFC matrix"""
    # Find perturbation column
    pert_col = None
    for col in ['condition', 'perturbation', 'perturb', 'pert_name']:
        if col in adata.obs.columns:
            pert_col = col
            break
    if pert_col is None:
        print(f"  [{cell_type_name}] No perturbation column found!")
        return None

    pert_values = list(adata.obs[pert_col].unique())

    # Identify control — exact match first
    ctrl_key = None
    for key in ['ctrl', 'control', 'NaN', 'nan', 'ctrl_ess']:
        if key in pert_values:
            ctrl_key = key
            break
    if ctrl_key is None:
        ctrl_key = adata.obs[pert_col].value_counts().index[0]

    print(f"  [{cell_type_name}] Control: '{ctrl_key}', pert_col: '{pert_col}', {len(pert_values)} unique values")

    # Compute control mean
    ctrl_mask = adata.obs[pert_col] == ctrl_key
    ctrl_mean = np.array(adata[ctrl_mask].X.mean(axis=0)).flatten()

    # Compute perturbation means
    perturbations = [p for p in pert_values if p != ctrl_key]
    actual_means = {}
    for pert in perturbations:
        mask = adata.obs[pert_col] == pert
        if mask.sum() > 0:
            actual_means[pert] = np.array(adata[mask].X.mean(axis=0)).flatten()

    # logFC matrix
    perts_with_data = sorted(actual_means.keys())
    true_logfc = np.array([actual_means[p] - ctrl_mean for p in perts_with_data])

    print(f"  [{cell_type_name}] {len(perts_with_data)} perturbations, {true_logfc.shape[1]} genes")
    print(f"  [{cell_type_name}] logFC: mean={true_logfc.mean():.4f}, std={true_logfc.std():.4f}, DEG_frac={((np.abs(true_logfc) > 0.25).mean()):.4f}")

    return {
        'cell_type': cell_type_name,
        'perturbations': perts_with_data,
        'true_logfc': true_logfc,
        'ctrl_mean': ctrl_mean,
        'n_perts': len(perts_with_data),
        'n_genes': true_logfc.shape[1],
    }

# ============================================================
# Model prediction generator
# ============================================================

def generate_model_predictions(true_logfc):
    """Generate 11 realistic model predictions"""
    gene_std = true_logfc.std(axis=0)
    median_std = np.median(gene_std[gene_std > 0])

    def gene_noise(scale, shape):
        base = np.random.normal(0, 1, shape)
        scale_vec = gene_std.reshape(1, -1) / (median_std + 1e-8)
        return base * scale_vec * scale

    np.random.seed(42)
    model_predictions = {}

    # 1. Mean predictor (baseline)
    model_predictions['mean_predictor'] = np.zeros_like(true_logfc)

    # 2. Additive linear (Ahlmann-Eltze best)
    model_predictions['additive_linear'] = true_logfc * 0.55 + gene_noise(0.03, true_logfc.shape)

    # 3. CPA-like
    model_predictions['cpa_like'] = true_logfc * 0.65 + gene_noise(0.06, true_logfc.shape)

    # 4. GEARS-like
    model_predictions['gears_like'] = true_logfc * 0.75 + gene_noise(0.08, true_logfc.shape)

    # 5. scGPT-like
    model_predictions['scgpt_like'] = true_logfc * 0.90 + gene_noise(0.04, true_logfc.shape)

    # 6. Well-calibrated noisy
    model_predictions['calibrated_noisy'] = true_logfc * 1.0 + gene_noise(0.12, true_logfc.shape)

    # 7. Over-predictor
    model_predictions['over_predictor'] = true_logfc * 1.3 + gene_noise(0.05, true_logfc.shape)

    # 8. Mean-effect trap
    mean_effect = true_logfc.mean(axis=1, keepdims=True)
    model_predictions['mean_effect_trap'] = mean_effect * 0.3 + gene_noise(0.02, true_logfc.shape)

    # 9. Partial sign-flip (20% of DEGs)
    deg_mask = np.abs(true_logfc) > 0.25
    partial_flip = true_logfc.copy()
    flip_mask = deg_mask & (np.random.random(true_logfc.shape) < 0.20)
    partial_flip[flip_mask] *= -1
    model_predictions['partial_flip'] = partial_flip * 0.8 + gene_noise(0.04, true_logfc.shape)

    # 10. Shuffled direction on DEGs
    shuffled = true_logfc.copy()
    deg_mask2 = np.abs(true_logfc) > 0.25
    random_sign = np.where(np.random.random(true_logfc.shape) > 0.5, 1.0, -1.0)
    shuffled[deg_mask2] = np.abs(shuffled[deg_mask2]) * random_sign[deg_mask2]
    model_predictions['shuffled_dir'] = shuffled * 0.6 + gene_noise(0.05, true_logfc.shape)

    # 11. Slightly-better-than-mean
    model_predictions['slight_above_mean'] = true_logfc * 0.15 + gene_noise(0.03, true_logfc.shape)

    return model_predictions

# ============================================================
# Main Analysis Pipeline
# ============================================================

def main():
    start_time = time.time()

    # ---- Phase 1: Load Data ----
    print("\n" + "=" * 60)
    print("PHASE 1: Data Loading + Model Predictions")
    print("=" * 60)

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
        print("ERROR: K562 data not found. Cannot proceed.")
        return

    # ---- Run analysis on each dataset ----
    all_dataset_results = {}

    for ct_name, ds_data in datasets.items():
        print(f"\n{'=' * 60}")
        print(f"ANALYZING: {ct_name}")
        print(f"{'=' * 60}")

        true_logfc = ds_data['true_logfc']
        n_perts = ds_data['n_perts']
        n_genes = ds_data['n_genes']

        # Generate model predictions
        model_predictions = generate_model_predictions(true_logfc)
        model_names = list(model_predictions.keys())
        print(f"Models ({len(model_names)}): {model_names}")

        # ---- Phase 2: Compute All Metrics ----
        print(f"\n  Phase 2: Computing metrics for {ct_name}...")
        deg_thresholds = [0.1, 0.25, 0.5, 1.0]
        all_metrics = {}

        for model_name in model_names:
            pred_mat = model_predictions[model_name]
            metrics = {}

            # Existing metrics
            existing = compute_existing_metrics(true_logfc, pred_mat)
            metrics['MSE'] = existing['mse_mean']
            metrics['R2'] = existing['r2_mean']
            metrics['Pearson'] = existing['pearson_mean']

            # BioEval-Dir (at multiple thresholds)
            for thresh in deg_thresholds:
                dir_result = bioeval_dir_from_matrices(true_logfc, pred_mat, deg_threshold=thresh)
                prefix = f"Dir_t{thresh}"
                metrics[f"{prefix}_all"] = dir_result['dir_accuracy_all']
                metrics[f"{prefix}_deg"] = dir_result['dir_accuracy_deg']
                metrics[f"{prefix}_weighted"] = dir_result['dir_accuracy_weighted']
                metrics[f"{prefix}_mag_ratio"] = dir_result['mag_ratio_all']
                metrics[f"{prefix}_mag_ratio_deg"] = dir_result['mag_ratio_deg']

            # BioEval-Cal
            cal_result = bioeval_cal(true_logfc, pred_mat)
            metrics['Cal_slope'] = cal_result['global_slope']
            metrics['Cal_intercept'] = cal_result['global_intercept']
            metrics['Cal_under_frac'] = cal_result['underpredict_frac']
            metrics['Cal_over_frac'] = cal_result['overpredict_frac']
            metrics['Cal_in_range'] = cal_result['cal_in_range_frac']

            # BioEval-DEG
            deg_result = bioeval_deg(true_logfc, pred_mat, deg_threshold=0.25)
            metrics['DEG_auprc'] = deg_result['auprc']
            metrics['DEG_dir_auprc'] = deg_result['dir_aware_auprc']
            metrics['DEG_jaccard'] = deg_result['jaccard']
            metrics['DEG_dir_jaccard'] = deg_result['dir_jaccard']

            all_metrics[model_name] = metrics

            # Print summary
            print(f"    {model_name}: MSE={metrics['MSE']:.4f}, R2={metrics['R2']:.4f}, "
                  f"Pearson={metrics['Pearson']:.3f}, Dir_all={metrics['Dir_t0.25_all']:.4f}, "
                  f"Dir_deg={metrics['Dir_t0.25_deg']:.4f}, Cal_slope={metrics['Cal_slope']:.3f}")

        # ---- Phase 3: Metric-Ranking Reversal Analysis ----
        print(f"\n  Phase 3: Ranking reversal analysis for {ct_name}...")

        ranking_metrics = {
            'MSE': {m: -all_metrics[m]['MSE'] for m in model_names},
            'R2': {m: all_metrics[m]['R2'] for m in model_names},
            'Pearson': {m: all_metrics[m]['Pearson'] for m in model_names},
            'Dir_all': {m: all_metrics[m]['Dir_t0.25_all'] for m in model_names},
            'Dir_deg': {m: all_metrics[m]['Dir_t0.25_deg'] for m in model_names},
            'Dir_weighted': {m: all_metrics[m]['Dir_t0.25_weighted'] for m in model_names},
            'Cal_slope_dev': {m: 1.0 - abs(all_metrics[m]['Cal_slope'] - 1.0) for m in model_names},
            'DEG_auprc': {m: all_metrics[m]['DEG_auprc'] for m in model_names},
            'DEG_dir_auprc': {m: all_metrics[m]['DEG_dir_auprc'] for m in model_names},
        }

        # Convert to rankings
        rankings = {}
        for metric_name, scores in ranking_metrics.items():
            sorted_models = sorted(scores.keys(), key=lambda m: scores[m], reverse=True)
            rankings[metric_name] = [sorted_models.index(m) + 1 for m in model_names]

        print(f"\n  Model Rankings ({ct_name}):")
        header = f"  {'Model':<22}"
        for mn in rankings:
            header += f"  {mn:>12}"
        print(header)
        for i, model in enumerate(model_names):
            row = f"  {model:<22}"
            for mn in rankings:
                row += f"  {rankings[mn][i]:>12}"
            print(row)

        # Kendall tau analysis — key pairs
        key_pairs = [
            ('MSE', 'Dir_all'),
            ('MSE', 'Dir_deg'),
            ('MSE', 'DEG_auprc'),
            ('MSE', 'DEG_dir_auprc'),
            ('R2', 'Dir_all'),
            ('R2', 'Dir_deg'),
            ('Pearson', 'Dir_all'),
            ('Pearson', 'Dir_deg'),
            ('Dir_all', 'Dir_deg'),
            ('DEG_auprc', 'DEG_dir_auprc'),
            ('Cal_slope_dev', 'Dir_all'),
            ('MSE', 'R2'),
        ]

        print(f"\n  Kendall tau ({ct_name}):")
        reversal_results = []
        for m1, m2 in key_pairs:
            tau, p = stats.kendalltau(rankings[m1], rankings[m2])
            interp = 'REVERSAL' if tau < 0.5 else ('PARTIAL' if tau < 0.7 else 'CONSISTENT')
            print(f"    tau({m1}, {m2}) = {tau:.3f} (p={p:.4f}) -> {interp}")
            reversal_results.append({
                'metric_1': m1, 'metric_2': m2,
                'kendall_tau': float(tau), 'p_value': float(p),
                'interpretation': interp
            })

        # Full pairwise
        full_kendall = kendall_rank_analysis(rankings)

        # ---- S1: DEG threshold sensitivity ----
        print(f"\n  S1: DEG threshold sensitivity for {ct_name}...")
        s1_results = {}
        for thresh in deg_thresholds:
            s1_results[f"t{thresh}"] = {}
            for model_name in model_names:
                dir_result = bioeval_dir_from_matrices(
                    true_logfc, model_predictions[model_name], deg_threshold=thresh
                )
                s1_results[f"t{thresh}"][model_name] = {
                    'dir_all': dir_result['dir_accuracy_all'],
                    'dir_deg': dir_result['dir_accuracy_deg'],
                    'mag_ratio': dir_result['mag_ratio_all'],
                }

        # Store dataset results
        all_dataset_results[ct_name] = {
            'data': {
                'cell_type': ct_name,
                'n_perturbations': n_perts,
                'n_genes': n_genes,
                'deg_threshold': 0.25,
                'frac_deg': float((np.abs(true_logfc) > 0.25).mean()),
            },
            'models': model_names,
            'all_metrics': all_metrics,
            'rankings': {k: v for k, v in rankings.items()},
            'reversal_analysis': reversal_results,
            'full_kendall': full_kendall,
            's1_threshold_sensitivity': s1_results,
        }

    # ---- Cross-cell-type comparison ----
    print(f"\n{'=' * 60}")
    print("CROSS-CELL-TYPE COMPARISON")
    print(f"{'=' * 60}")

    cross_ct_results = {}
    if 'K562' in all_dataset_results and 'RPE1' in all_dataset_results:
        k562_rank = all_dataset_results['K562']['rankings']
        rpe1_rank = all_dataset_results['RPE1']['rankings']
        model_names_k562 = all_dataset_results['K562']['models']
        model_names_rpe1 = all_dataset_results['RPE1']['models']

        # Compare rankings across cell types for shared metrics
        shared_metrics = set(k562_rank.keys()) & set(rpe1_rank.keys())
        print(f"\nShared metrics between K562 and RPE1: {shared_metrics}")

        # Only compare if same models
        if model_names_k562 == model_names_rpe1:
            print("\nSame models -- comparing cross-cell-type ranking consistency:")
            for metric in sorted(shared_metrics):
                tau, p = stats.kendalltau(k562_rank[metric], rpe1_rank[metric])
                interp = 'CONSISTENT' if tau > 0.7 else ('PARTIAL' if tau > 0.5 else 'VARIABLE')
                print(f"  {metric}: tau(K562, RPE1) = {tau:.3f} (p={p:.4f}) -> {interp}")
                cross_ct_results[metric] = {'tau': float(tau), 'p': float(p), 'interpretation': interp}
        else:
            print(f"  Different models: K562 has {len(model_names_k562)}, RPE1 has {len(model_names_rpe1)}")

    # ---- Compile Final Results ----
    elapsed = time.time() - start_time

    final_results = {
        'run': 'run_13_v2',
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'runtime_s': elapsed,
        'n_models': len(list(all_dataset_results.values())[0]['models']) if all_dataset_results else 0,
        'datasets': all_dataset_results,
        'cross_cell_type': cross_ct_results,
    }

    results_path = os.path.join(OUTPUT_DIR, "run_13_results.json")
    with open(results_path, 'w') as f:
        json.dump(final_results, f, indent=2, default=str)
    print(f"\nResults saved to: {results_path}")

    # ---- Summary ----
    print(f"\n{'=' * 60}")
    print("SUMMARY")
    print(f"{'=' * 60}")
    for ct_name, ct_results in all_dataset_results.items():
        reversals = [r for r in ct_results['reversal_analysis'] if r['interpretation'] == 'REVERSAL']
        partials = [r for r in ct_results['reversal_analysis'] if r['interpretation'] == 'PARTIAL']
        consistent = [r for r in ct_results['reversal_analysis'] if r['interpretation'] == 'CONSISTENT']
        print(f"\n  {ct_name} ({ct_results['data']['n_perturbations']} perts, {ct_results['data']['n_genes']} genes):")
        print(f"    REVERSAL: {len(reversals)}, PARTIAL: {len(partials)}, CONSISTENT: {len(consistent)}")
        if reversals:
            print(f"    Key reversals:")
            for r in reversals:
                print(f"      tau({r['metric_1']}, {r['metric_2']}) = {r['kendall_tau']:.3f}")

    print(f"\nTotal runtime: {elapsed:.1f}s")
    print("Analysis complete!")

if __name__ == "__main__":
    main()
