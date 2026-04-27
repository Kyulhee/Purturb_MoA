"""
NAP Transfer Experiment v2 — Main Entry Point
================================================
Runs the full NAP transfer experiment with common embedding space.
"""

import time
import warnings
import numpy as np
from typing import Dict

import cobra
import xgboost as xgb
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score

warnings.filterwarnings("ignore")

# Import from the companion module
from nap_transfer_experiment_v2 import (
    load_model, generate_fba_dataset, model_to_heterodata,
    HGTGNNEncoder, SharedProjectionHead, GrowthRegressionHead,
    MetabolicSurrogateModel,
    train_end_to_end, train_encoder_only,
    extract_common_embeddings,
    experiment_xgboost_only, experiment_gnn_end_to_end,
    experiment_transfer, experiment_gnn_xgboost_hybrid,
    analyze_feature_completeness, analyze_embedding_alignment,
)


def run_nap_experiment():
    print("=" * 70)
    print("NAP Transfer Experiment v2: Common Embedding Space")
    print("=" * 70)

    # --- Step 1: Load models and generate data ---
    print("\n[Step 1] Loading models and generating FBA data...")

    models_config = [
        ("textbook", "E. coli core", 137),
        ("iSB619", "S. aureus", 619),
        ("iCN718", "C. neoformans", 709),
    ]

    datasets = {}
    graphs = {}

    for model_id, organism, n_genes in models_config:
        print(f"\n  --- {model_id} ({organism}) ---")
        model = load_model(model_id)
        wt = model.optimize().objective_value
        print(f"  WT growth: {wt:.4f}")

        n_samples = 300 if model_id == "textbook" else 150
        masks, growth = generate_fba_dataset(model, n_samples=n_samples, seed=42)
        datasets[model_id] = {"masks": masks, "growth": growth, "n_genes": n_genes}

        graph = model_to_heterodata(model)
        graphs[model_id] = graph
        n_nodes = (graph["metabolite"].num_nodes +
                   graph["reaction"].num_nodes +
                   graph["gene"].num_nodes)
        print(f"  Graph: {n_nodes} nodes ({graph['metabolite'].num_nodes} met + "
              f"{graph['reaction'].num_nodes} rxn + {graph['gene'].num_nodes} gene)")

    # --- Step 2: XGBoost-only baseline ---
    print("\n[Step 2] XGBoost-only baselines (5-fold CV)...")

    xgb_results = {}
    for model_id in datasets:
        masks = datasets[model_id]["masks"]
        growth = datasets[model_id]["growth"]
        r2_mean, r2_std = experiment_xgboost_only(masks, growth, n_folds=5)
        xgb_results[model_id] = {"r2_mean": r2_mean, "r2_std": r2_std}
        print(f"  {model_id}: R2 = {r2_mean:.4f} +/- {r2_std:.4f}")

    # --- Step 3: GNN End-to-End on each model ---
    print("\n[Step 3] GNN End-to-End (SharedProjection + RegHead, 3-fold CV)...")

    e2e_results = {}
    trained_models = {}

    for model_id in datasets:
        print(f"\n  --- {model_id} ---")
        masks = datasets[model_id]["masks"]
        growth = datasets[model_id]["growth"]
        graph = graphs[model_id]

        try:
            r2_mean, r2_std = experiment_gnn_end_to_end(
                graph, masks, growth, n_folds=3,
                hidden=32, out=32, common_dim=16,
                epochs=30, lr=0.005
            )
            e2e_results[model_id] = {"r2_mean": r2_mean, "r2_std": r2_std}
            print(f"  {model_id}: E2E R2 = {r2_mean:.4f} +/- {r2_std:.4f}")

            # Train a full model on ALL data for transfer
            metadata = graph.metadata()
            encoder = HGTGNNEncoder(metadata, hidden_channels=32,
                                     out_channels=32, num_heads=2, num_layers=2)
            projection = SharedProjectionHead(gnn_dim=32, common_dim=16)
            regression = GrowthRegressionHead(common_dim=16)
            full_model = MetabolicSurrogateModel(encoder, projection, regression)
            full_model = train_end_to_end(
                full_model, graph, masks, growth, epochs=30, lr=0.005, verbose=True)
            trained_models[model_id] = full_model

        except Exception as e:
            print(f"  {model_id}: FAILED - {e}")
            import traceback
            traceback.print_exc()
            e2e_results[model_id] = {"r2_mean": -999, "r2_std": 0}

    # --- Step 4: Cross-species Transfer via Common Embedding Space ---
    print("\n[Step 4] Cross-species Transfer (Common Embedding Space)...")

    transfer_results = {}
    transfer_pairs = [
        ("textbook", "iSB619"),
        ("textbook", "iCN718"),
    ]

    for source_id, target_id in transfer_pairs:
        if source_id not in trained_models:
            print(f"  SKIP {source_id}->{target_id}: source model not available")
            continue

        print(f"\n  --- {source_id} -> {target_id} ---")
        source_model = trained_models[source_id]
        target_masks = datasets[target_id]["masks"]
        target_growth = datasets[target_id]["growth"]
        target_graph = graphs[target_id]

        try:
            r2_scratch_mean, r2_scratch_std, r2_transfer_mean, r2_transfer_std = \
                experiment_transfer(
                    source_model, target_graph, target_masks, target_growth,
                    n_folds=3, fine_tune_epochs=30, lr=0.005, common_dim=16
                )
            transfer_results[f"{source_id}->{target_id}"] = {
                "scratch_r2_mean": r2_scratch_mean,
                "scratch_r2_std": r2_scratch_std,
                "transfer_r2_mean": r2_transfer_mean,
                "transfer_r2_std": r2_transfer_std,
                "delta": r2_transfer_mean - r2_scratch_mean,
            }
            print(f"  Scratch R2 = {r2_scratch_mean:.4f} +/- {r2_scratch_std:.4f}")
            print(f"  Transfer R2 = {r2_transfer_mean:.4f} +/- {r2_transfer_std:.4f}")
            print(f"  Delta = {r2_transfer_mean - r2_scratch_mean:+.4f}")
        except Exception as e:
            print(f"  FAILED: {e}")
            import traceback
            traceback.print_exc()

    # --- Step 5: Embedding Alignment Analysis (C3 validation) ---
    print("\n[Step 5] Common Embedding Space Alignment (C3 condition)...")

    for source_id, target_id in transfer_pairs:
        if source_id not in trained_models or target_id not in trained_models:
            continue
        print(f"\n  {source_id} <-> {target_id}:")
        try:
            mean_corr, per_comp = analyze_embedding_alignment(
                trained_models[source_id], trained_models[target_id],
                graphs[source_id], graphs[target_id],
                datasets[source_id]["masks"], datasets[target_id]["masks"],
                n_samples=50
            )
        except Exception as e:
            print(f"  Alignment analysis failed: {e}")

    # --- Step 6: Feature Completeness Analysis (C1 condition) ---
    print("\n[Step 6] Feature completeness analysis (C1 condition)...")

    for model_id in datasets:
        print(f"\n  {model_id}:")
        masks = datasets[model_id]["masks"]
        growth = datasets[model_id]["growth"]
        ratio = analyze_feature_completeness(masks, growth)

    # --- Step 7: Graph Structural Comparison (C2 condition) ---
    print("\n[Step 7] Graph structural comparison (C2 condition)...")

    model_ids = list(graphs.keys())
    for i in range(len(model_ids)):
        for j in range(i + 1, len(model_ids)):
            g1, g2 = graphs[model_ids[i]], graphs[model_ids[j]]
            n1 = {nt: g1[nt].num_nodes for nt in ["metabolite", "reaction", "gene"]}
            n2 = {nt: g2[nt].num_nodes for nt in ["metabolite", "reaction", "gene"]}
            ratio = {nt: n2[nt] / max(n1[nt], 1) for nt in n1}
            different = any(r != 1.0 for r in ratio.values())
            print(f"  {model_ids[i]} vs {model_ids[j]}: "
                  f"met {n1['metabolite']}->{n2['metabolite']}, "
                  f"rxn {n1['reaction']}->{n2['reaction']}, "
                  f"gene {n1['gene']}->{n2['gene']} | "
                  f"C2(T)={'T' if different else 'F'}")

    # --- Step 8: Summary ---
    print("\n" + "=" * 70)
    print("NAP PREDICTION vs EXPERIMENTAL RESULT")
    print("=" * 70)

    print(f"\n{'Model':<15} {'NAP':<8} {'Pred':<18} "
          f"{'XGB R2':<12} {'GNN E2E R2':<12} {'GNN>XGB?':<10}")
    print("-" * 75)

    nap_scores = {
        "textbook": 0,  # C1=F, C2=F (fixed graph within species)
        "iSB619": 0,    # Same: within single species, C2=F
        "iCN718": 0,    # Same
    }

    for model_id in datasets:
        xgb_r2 = xgb_results.get(model_id, {}).get("r2_mean", -999)
        gnn_r2 = e2e_results.get(model_id, {}).get("r2_mean", -999)
        nap = nap_scores.get(model_id, 0)
        pred = "GNN value" if nap >= 2 else "GNN no value"
        gnn_better = "YES" if gnn_r2 > xgb_r2 + 0.02 else (
            "MARGINAL" if gnn_r2 > xgb_r2 else "NO")
        print(f"{model_id:<15} {nap}/6{'':<4} {pred:<18} "
              f"{xgb_r2:<12.4f} {gnn_r2:<12.4f} {gnn_better:<10}")

    print(f"\n{'Transfer':<25} {'NAP':<8} {'Pred':<18} "
          f"{'Scratch R2':<12} {'Transfer R2':<12} {'Delta':<10}")
    print("-" * 85)

    for key, res in transfer_results.items():
        source_id, target_id = key.split("->")
        nap = 2  # C2=T + C3=T for cross-species transfer
        pred = "GNN value expected"
        print(f"{key:<25} {nap}/6{'':<4} {pred:<18} "
              f"{res['scratch_r2_mean']:<12.4f} {res['transfer_r2_mean']:<12.4f} "
              f"{res['delta']:+.4f}")

    # --- NAP Validation Summary ---
    print("\n" + "=" * 70)
    print("NAP VALIDATION SUMMARY")
    print("=" * 70)

    print("\n1. C1 (Feature Incompleteness):")
    for model_id in datasets:
        masks = datasets[model_id]["masks"]
        growth = datasets[model_id]["growth"]
        from sklearn.feature_selection import mutual_info_regression
        mi = mutual_info_regression(masks, growth, random_state=42)
        ratio = mi.sum() / np.var(growth) if np.var(growth) > 0 else 0
        c1 = "T" if ratio < 0.9 else "F"
        print(f"  {model_id}: C1={c1} (MI/Var={ratio:.4f})")

    print("\n2. C2 (Graph Variability):")
    print("  Within single species: C2=F (graph is fixed)")
    print("  Across species: C2=T (different graph structures)")

    print("\n3. C3 (Transfer Learning):")
    for key, res in transfer_results.items():
        delta = res["delta"]
        c3_valid = "T" if delta > 0.02 else "F"
        print(f"  {key}: C3={c3_valid} (transfer delta={delta:+.4f})")

    print("\n4. Key Insight:")
    print("  The common embedding space approach makes cross-species transfer")
    print("  possible by decoupling species-specific GNN encoders from the")
    print("  shared projection + regression head.")
    print("  Transfer success depends on whether the learned common space")
    print("  captures universal metabolic principles that generalize across species.")

    return {
        "xgb_results": xgb_results,
        "e2e_results": e2e_results,
        "transfer_results": transfer_results,
        "datasets": {k: {"n_genes": v["n_genes"],
                         "n_samples": len(v["growth"])}
                     for k, v in datasets.items()},
    }


if __name__ == "__main__":
    results = run_nap_experiment()
