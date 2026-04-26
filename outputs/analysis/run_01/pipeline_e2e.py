"""
End-to-End Pipeline: Module A → B → C
=======================================
Full pipeline for GNN surrogate model + Active Learning
for FBA-based metabolic network optimization.

Research Question:
  "대사 네트워크의 이종 그래프에서 GNN 대리 모델이 FBA 기반 탐색 공간 축소에 기여하는가?"

Evaluation:
  1. GNN embedding effect: XGBoost-only vs GNN+XGBoost (R2 comparison)
  2. AL efficiency: Random vs AL (FBA calls to reach same R2)
  3. AL transition: diversity-only vs two-phase (R2 convergence curve)
  4. GNN pretraining: no-pretrain vs edge-prediction (R2)
"""

import sys
import os
import time
import json
import numpy as np
import torch
from datetime import datetime

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from module_a_fba_generator import (
    FBAGroundTruthGenerator,
    generate_single_knockouts,
    generate_double_knockouts,
    generate_random_knockouts,
    build_knockout_mask,
    _run_knockout_batch,
    model_to_heterodata,
)
from module_b_gnn_xgboost_surrogate import (
    GNNXGBoostSurrogate,
    XGBoostOnlyBaseline,
    HGTGNN,
)
from module_c_active_learning import (
    ActiveLearningLoop,
    RandomScreeningBaseline,
    generate_candidate_pool,
)


def fba_oracle(masks: np.ndarray, model, n_workers: int = 1) -> np.ndarray:
    """
    FBA oracle: given knockout masks, return growth rates.
    Wraps Module A's FBA execution.
    """
    gene_ids = sorted([g.id for g in model.genes])
    gene_idx = {gid: i for i, gid in enumerate(gene_ids)}

    # Convert masks to gene ID lists
    knockout_combos = []
    for i in range(len(masks)):
        ko_indices = np.where(masks[i] > 0.5)[0]
        ko_genes = [gene_ids[j] for j in ko_indices if j < len(gene_ids)]
        knockout_combos.append(ko_genes)

    growth_rates = _run_knockout_batch(model, knockout_combos, n_workers=n_workers)
    return np.array(growth_rates, dtype=np.float32)


def run_experiment_1_gnn_effect(data, model, verbose=True):
    """
    Experiment 1: GNN embedding effect.
    Compare XGBoost-only vs GNN+XGBoost.
    """
    print("\n" + "=" * 60)
    print("EXPERIMENT 1: GNN Embedding Effect")
    print("=" * 60)

    # Combine all knockout data
    all_masks = []
    all_growth = []

    for key_mask, key_growth in [
        ("single_ko_mask", "single_ko_growth"),
        ("double_ko_mask", "double_ko_growth"),
        ("random_ko_mask", "random_ko_growth"),
    ]:
        if key_mask in data:
            mask = data[key_mask]
            if isinstance(mask, torch.Tensor):
                mask = mask.numpy()
            all_masks.append(mask)
            all_growth.append(data[key_growth])

    X = np.concatenate(all_masks, axis=0)
    y = np.concatenate(all_growth, axis=0)

    if verbose:
        print(f"  Dataset: {X.shape[0]} samples, {X.shape[1]} features")
        print(f"  Growth: [{y.min():.4f}, {y.max():.4f}], mean={y.mean():.4f}")

    results = {}

    # XGBoost-only baseline
    if verbose:
        print("\n--- XGBoost-only Baseline ---")
    baseline = XGBoostOnlyBaseline()
    baseline_results = baseline.fit(X, y, verbose=verbose)
    results["xgb_only"] = baseline_results

    # GNN + XGBoost
    if verbose:
        print("\n--- GNN + XGBoost ---")
    graph = data["graph"]
    surrogate = GNNXGBoostSurrogate(
        graph,
        hidden_channels=32,
        out_channels=32,
        num_heads=2,
        num_layers=2,
    )
    gnn_results = surrogate.fit(X, y, pretrain_epochs=30, verbose=verbose)
    results["gnn_xgb"] = gnn_results

    # Compare
    if verbose:
        print(f"\n  R2 comparison:")
        print(f"    XGBoost-only: {baseline_results['r2_test']:.4f}")
        print(f"    GNN+XGBoost:  {gnn_results['r2_test']:.4f}")
        improvement = gnn_results['r2_test'] - baseline_results['r2_test']
        print(f"    Improvement:  {improvement:+.4f}")

    return results


def run_experiment_2_al_efficiency(data, model, verbose=True):
    """
    Experiment 2: AL efficiency.
    Compare Random screening vs Active Learning (FBA calls to reach same R2).
    """
    print("\n" + "=" * 60)
    print("EXPERIMENT 2: AL Efficiency")
    print("=" * 60)

    # Use single + small random as initial data
    initial_masks = data["single_ko_mask"]
    if isinstance(initial_masks, torch.Tensor):
        initial_masks = initial_masks.numpy()
    initial_growth = data["single_ko_growth"]

    n_genes = initial_masks.shape[1]

    def fba_fn(masks):
        return fba_oracle(masks, model, n_workers=1)

    def surrogate_factory():
        return GNNXGBoostSurrogate(
            data["graph"],
            hidden_channels=16,
            out_channels=16,
            num_heads=2,
            num_layers=2,
        )

    results = {}

    # Active Learning
    if verbose:
        print("\n--- Active Learning ---")
    al = ActiveLearningLoop(
        n_genes=n_genes,
        n_select_per_round=50,
        max_rounds=10,  # Reduced for testing
        transition_r2=0.3,
        transition_patience=3,
        ucb_alpha=0.5,
    )
    al_results = al.run(
        initial_masks=initial_masks,
        initial_growth=initial_growth,
        fba_oracle=fba_fn,
        surrogate_factory=surrogate_factory,
        verbose=verbose,
    )
    results["al"] = {
        "total_fba_calls": al_results["total_fba_calls"],
        "final_r2": al_results["final_r2"],
        "total_rounds": al_results["total_rounds"],
        "phase_transition_round": al_results["phase_transition_round"],
        "r2_per_round": al_results["r2_per_round"],
    }

    # Random Screening (same budget)
    if verbose:
        print("\n--- Random Screening ---")
    random_bl = RandomScreeningBaseline(
        n_genes=n_genes,
        n_select_per_round=50,
        max_rounds=10,
    )
    random_results = random_bl.run(
        initial_masks=initial_masks,
        initial_growth=initial_growth,
        fba_oracle=fba_fn,
        surrogate_factory=surrogate_factory,
        verbose=verbose,
    )
    results["random"] = {
        "total_fba_calls": random_results["total_fba_calls"],
        "final_r2": random_results["final_r2"],
        "total_rounds": random_results["total_rounds"],
        "r2_per_round": random_results["r2_per_round"],
    }

    # Compare
    if verbose:
        print(f"\n  FBA call comparison:")
        print(f"    AL:     {results['al']['total_fba_calls']} calls, R2={results['al']['final_r2']:.4f}")
        print(f"    Random: {results['random']['total_fba_calls']} calls, R2={results['random']['final_r2']:.4f}")

    return results


def run_experiment_3_pretraining(data, model, verbose=True):
    """
    Experiment 3: GNN pretraining comparison.
    No-pretrain vs edge-prediction pretraining.
    """
    print("\n" + "=" * 60)
    print("EXPERIMENT 3: GNN Pretraining Comparison")
    print("=" * 60)

    all_masks = []
    all_growth = []
    for key_mask, key_growth in [
        ("single_ko_mask", "single_ko_growth"),
        ("double_ko_mask", "double_ko_growth"),
        ("random_ko_mask", "random_ko_growth"),
    ]:
        if key_mask in data:
            mask = data[key_mask]
            if isinstance(mask, torch.Tensor):
                mask = mask.numpy()
            all_masks.append(mask)
            all_growth.append(data[key_growth])

    X = np.concatenate(all_masks, axis=0)
    y = np.concatenate(all_growth, axis=0)

    graph = data["graph"]
    results = {}

    # No pretraining
    if verbose:
        print("\n--- No Pretraining ---")
    surrogate_no_pt = GNNXGBoostSurrogate(
        graph, hidden_channels=16, out_channels=16
    )
    no_pt_results = surrogate_no_pt.fit(X, y, pretrain_epochs=0, verbose=verbose)
    results["no_pretrain"] = no_pt_results

    # Edge prediction pretraining
    if verbose:
        print("\n--- Edge Prediction Pretraining ---")
    surrogate_edge_pt = GNNXGBoostSurrogate(
        graph, hidden_channels=16, out_channels=16
    )
    edge_pt_results = surrogate_edge_pt.fit(X, y, pretrain_epochs=30, verbose=verbose)
    results["edge_pretrain"] = edge_pt_results

    if verbose:
        print(f"\n  R2 comparison:")
        print(f"    No pretrain:     {no_pt_results['r2_test']:.4f}")
        print(f"    Edge prediction: {edge_pt_results['r2_test']:.4f}")

    return results


def run_full_pipeline(verbose=True):
    """Run the complete end-to-end pipeline."""

    print("=" * 60)
    print("END-TO-END PIPELINE: Module A -> B -> C")
    print(f"Started: {datetime.now().isoformat()}")
    print("=" * 60)

    total_t0 = time.time()
    all_results = {}

    # ── Module A: Generate FBA ground truth ──
    print("\n" + "=" * 60)
    print("MODULE A: FBA Ground Truth Generation")
    print("=" * 60)

    t0 = time.time()
    gen = FBAGroundTruthGenerator()
    data = gen.run(
        single=True,
        double=True,
        double_n=500,
        random_n=200,
        n_workers=1,
    )
    module_a_time = time.time() - t0
    all_results["module_a_time"] = module_a_time
    all_results["wt_growth"] = data["wt_growth"]
    all_results["n_samples"] = sum(
        len(data[k]) for k in ["single_ko_growth", "double_ko_growth", "random_ko_growth"]
        if k in data
    )

    if verbose:
        print(f"  Total samples: {all_results['n_samples']}")
        print(f"  Time: {module_a_time:.1f}s")

    # ── Experiments ──
    try:
        exp1_results = run_experiment_1_gnn_effect(data, gen.model, verbose=verbose)
        all_results["exp1_gnn_effect"] = {
            "xgb_only_r2": exp1_results["xgb_only"]["r2_test"],
            "gnn_xgb_r2": exp1_results["gnn_xgb"]["r2_test"],
        }
    except Exception as e:
        print(f"  Experiment 1 FAILED: {e}")
        all_results["exp1_gnn_effect"] = {"error": str(e)}

    try:
        exp2_results = run_experiment_2_al_efficiency(data, gen.model, verbose=verbose)
        all_results["exp2_al_efficiency"] = {
            "al_fba_calls": exp2_results["al"]["total_fba_calls"],
            "al_final_r2": exp2_results["al"]["final_r2"],
            "random_fba_calls": exp2_results["random"]["total_fba_calls"],
            "random_final_r2": exp2_results["random"]["final_r2"],
        }
    except Exception as e:
        print(f"  Experiment 2 FAILED: {e}")
        all_results["exp2_al_efficiency"] = {"error": str(e)}

    try:
        exp3_results = run_experiment_3_pretraining(data, gen.model, verbose=verbose)
        all_results["exp3_pretraining"] = {
            "no_pretrain_r2": exp3_results["no_pretrain"]["r2_test"],
            "edge_pretrain_r2": exp3_results["edge_pretrain"]["r2_test"],
        }
    except Exception as e:
        print(f"  Experiment 3 FAILED: {e}")
        all_results["exp3_pretraining"] = {"error": str(e)}

    # ── Summary ──
    total_time = time.time() - total_t0
    all_results["total_time"] = total_time

    print("\n" + "=" * 60)
    print("PIPELINE SUMMARY")
    print("=" * 60)
    print(f"  Total time: {total_time:.1f}s")
    print(f"  Module A time: {module_a_time:.1f}s")
    print(f"  WT growth: {data['wt_growth']:.4f}")
    print(f"  Total FBA samples: {all_results['n_samples']}")

    if "exp1_gnn_effect" in all_results and "error" not in all_results["exp1_gnn_effect"]:
        r = all_results["exp1_gnn_effect"]
        print(f"\n  Exp1 - GNN Effect:")
        print(f"    XGBoost-only R2: {r['xgb_only_r2']:.4f}")
        print(f"    GNN+XGBoost R2:  {r['gnn_xgb_r2']:.4f}")

    if "exp2_al_efficiency" in all_results and "error" not in all_results["exp2_al_efficiency"]:
        r = all_results["exp2_al_efficiency"]
        print(f"\n  Exp2 - AL Efficiency:")
        print(f"    AL FBA calls: {r['al_fba_calls']}, R2: {r['al_final_r2']:.4f}")
        print(f"    Random FBA calls: {r['random_fba_calls']}, R2: {r['random_final_r2']:.4f}")

    if "exp3_pretraining" in all_results and "error" not in all_results["exp3_pretraining"]:
        r = all_results["exp3_pretraining"]
        print(f"\n  Exp3 - Pretraining:")
        print(f"    No pretrain R2:     {r['no_pretrain_r2']:.4f}")
        print(f"    Edge prediction R2: {r['edge_pretrain_r2']:.4f}")

    print(f"\n  Completed: {datetime.now().isoformat()}")

    # Save results
    results_path = os.path.join(os.path.dirname(__file__), "pipeline_results.json")
    # Convert numpy types to Python types for JSON
    json_results = {}
    for k, v in all_results.items():
        if isinstance(v, (np.floating, np.integer)):
            json_results[k] = float(v)
        elif isinstance(v, dict):
            json_results[k] = {
                kk: float(vv) if isinstance(vv, (np.floating, np.integer)) else vv
                for kk, vv in v.items()
            }
        else:
            json_results[k] = v

    with open(results_path, "w", encoding="utf-8") as f:
        json.dump(json_results, f, indent=2, ensure_ascii=False, default=str)
    print(f"\n  Results saved to: {results_path}")

    return all_results


if __name__ == "__main__":
    run_full_pipeline()
