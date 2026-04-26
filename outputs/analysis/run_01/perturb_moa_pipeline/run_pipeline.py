"""
Cross-Perturbation MoA Framework v2 -- Main Runner

Key changes from v1:
- Uses rFCFP fingerprint-based drug encoding
- Evaluates metric (prototype) head + zero-shot nearest-neighbor
- Compares classification head vs metric head vs zero-shot NN

Usage:
    python run_pipeline.py --mode [full|quick|eval]
"""
import sys
import argparse
import numpy as np
import torch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from config import ModelConfig, RESULTS_DIR, SEED
from data_utils import (
    create_synthetic_sci_plex, annotate_moa, preprocess_sci_plex,
    prepare_dataloaders
)
from model import CrossPerturbMoAModel
from train import train_pipeline, evaluate_cross_perturbation, set_seed


def run_synthetic_quick():
    """Quick pipeline test with synthetic data + rFCFP encoding."""
    print("=" * 80)
    print("Cross-Perturbation MoA Framework v2 -- Synthetic Quick Test")
    print("  - rFCFP fingerprint-based drug encoding")
    print("  - Metric (prototype) head for zero-shot MoA")
    print("  - Nearest-neighbor zero-shot evaluation")
    print("=" * 80)

    set_seed(SEED)

    # 1. Create synthetic data
    print("\n[1/4] Creating synthetic sci-Plex-like data...")
    adata = create_synthetic_sci_plex(
        n_cells=20000,
        n_genes=1000,
        n_drugs=30,
        n_cell_lines=3,
        n_moa_classes=6,
    )
    adata = annotate_moa(adata)

    # 2. Prepare data splits
    print("\n[2/4] Preparing data splits with fingerprints...")
    strategies = ["random", "leave_compound_out", "leave_moa_out"]
    results_all = {}

    for strategy in strategies:
        print(f"\n--- Split strategy: {strategy} ---")
        try:
            train_ds, test_ds, meta = prepare_dataloaders(
                adata, split_strategy=strategy, seed=SEED,
                fingerprint_bits=512,  # smaller for quick test
            )
        except Exception as e:
            print(f"  Skipped ({e})")
            import traceback
            traceback.print_exc()
            continue

        if meta['train_size'] < 100 or meta['test_size'] < 50:
            print(f"  Skipped (too few samples: train={meta['train_size']}, test={meta['test_size']})")
            continue

        # Adjust config for synthetic data
        config = ModelConfig()
        config.n_epochs = 30
        config.eval_every = 5
        config.batch_size = 128
        config.learning_rate = 1e-3
        config.latent_dim = 64
        config.basal_dim = 32
        config.perturbation_dim = 16
        config.covariate_dim = 8
        config.drug_hidden_dims = [128, 64]  # smaller MLP for quick test
        config.fingerprint_bits = 512

        # 3. Train
        save_dir = RESULTS_DIR / f"v2_synthetic_{strategy}"
        print(f"\n[3/4] Training ({strategy})...")
        try:
            model, results, z_test, moa_true = train_pipeline(
                train_ds, test_ds, meta, config, save_dir
            )
        except Exception as e:
            print(f"  Training failed: {e}")
            import traceback
            traceback.print_exc()
            continue

        # 4. Cross-perturbation evaluation
        print(f"\n[4/4] Cross-perturbation matching...")
        try:
            cp_results = evaluate_cross_perturbation(
                model, drug_z=z_test, drug_moa=moa_true, top_k=5
            )
            results.update(cp_results)
        except Exception as e:
            print(f"  Cross-perturbation eval failed: {e}")

        results_all[strategy] = results

        # Print summary
        print(f"\n{'='*70}")
        print(f"Strategy: {strategy}")
        print(f"{'='*70}")
        print(f"  Cls Head Top-1:      {results.get('cls_accuracy_top1', 0):.4f}")
        print(f"  Metric Head Top-1:   {results.get('metric_accuracy_top1', 0):.4f}")
        print(f"  Zero-shot NN Top-1:  {results.get('zero_shot_top1_acc', 0):.4f}")
        print(f"  Recon MSE:           {results.get('recon_mse', 0):.4f}")
        print(f"  Cross-Perturb:       {results.get('within_drug_top1_match', 0):.4f}")

    # Print overall summary
    print(f"\n{'='*80}")
    print("SUMMARY -- v2 Synthetic Quick Test (rFCFP + Metric Head)")
    print(f"{'='*80}")
    print(f"{'Strategy':25s} | {'Cls Top-1':>10s} | {'Metric Top-1':>12s} | "
          f"{'ZeroShot Top-1':>14s} | {'Recon MSE':>10s}")
    print("-" * 80)
    for strategy, res in results_all.items():
        print(f"  {strategy:23s} | {res.get('cls_accuracy_top1', 0):10.4f} | "
              f"{res.get('metric_accuracy_top1', 0):12.4f} | "
              f"{res.get('zero_shot_top1_acc', 0):14.4f} | "
              f"{res.get('recon_mse', 0):10.4f}")

    return results_all


def run_full_pipeline():
    """Full pipeline with real sci-Plex data."""
    print("=" * 80)
    print("Cross-Perturbation MoA Framework v2 -- Full Pipeline")
    print("=" * 80)

    print("\nReal data pipeline requires sci-Plex download.")
    print("Falling back to synthetic quick test...")
    return run_synthetic_quick()


def main():
    parser = argparse.ArgumentParser(description="Cross-Perturbation MoA Framework v2")
    parser.add_argument("--mode", choices=["full", "quick", "eval"],
                        default="quick", help="Pipeline mode")
    parser.add_argument("--seed", type=int, default=SEED, help="Random seed")
    args = parser.parse_args()

    set_seed(args.seed)

    if args.mode == "quick":
        results = run_synthetic_quick()
    elif args.mode == "full":
        results = run_full_pipeline()
    elif args.mode == "eval":
        print("Evaluation mode: load model from results/best_model.pt")

    print("\nDone.")


if __name__ == "__main__":
    main()
