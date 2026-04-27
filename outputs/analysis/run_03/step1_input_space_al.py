"""
Step 1: Input-space Active Learning Micro-Validation
=====================================================
Goal: Verify that input-space AL (Hamming distance + ensemble uncertainty)
outperforms random sampling for XGBoost surrogate.

Strategy:
  - Phase 1 (R2 < 0.3): diversity-based -- max Hamming distance from selected
  - Phase 2 (R2 >= 0.3): UCB -- ensemble uncertainty + predicted growth

Micro-validation: 200 initial -> 5 rounds x 20 samples = 100 AL samples
Compare against: 200 initial -> 5 rounds x 20 random samples
"""

import copy
import itertools
import random
import time
from typing import Dict, List, Optional, Tuple

import cobra
import numpy as np
import xgboost as xgb
from sklearn.model_selection import cross_val_score
from scipy.spatial.distance import cdist


# ============================================================
# 1. FBA Data Generation
# ============================================================

def load_model():
    model = cobra.io.load_model("textbook")
    return model

def run_fba_knockout(model, gene_ids):
    """Single FBA knockout, return growth rate."""
    try:
        model_cp = model.copy()
        with model_cp:
            for gid in gene_ids:
                if gid in [g.id for g in model_cp.genes]:
                    model_cp.genes.get_by_id(gid).knock_out()
            sol = model_cp.optimize()
            if sol.status == "optimal":
                return sol.objective_value
        return 0.0
    except:
        return 0.0

def build_knockout_mask(model, combos):
    """Build binary mask (n_combos, n_genes)."""
    gene_ids = sorted([g.id for g in model.genes])
    gene_idx = {gid: i for i, gid in enumerate(gene_ids)}
    n_gene = len(gene_ids)
    mask = np.zeros((len(combos), n_gene), dtype=np.float32)
    for i, combo in enumerate(combos):
        for gid in combo:
            if gid in gene_idx:
                mask[i, gene_idx[gid]] = 1.0
    return mask, gene_ids

def generate_random_combos(model, n, min_k=1, max_k=5, seed=42):
    """Generate random knockout combinations."""
    gene_ids = [g.id for g in model.genes]
    rng = random.Random(seed)
    combos = []
    for _ in range(n):
        k = rng.randint(min_k, min(max_k, len(gene_ids)))
        combo = rng.sample(gene_ids, k)
        combos.append(combo)
    return combos

def compute_growth_batch(model, combos):
    """Compute growth rates for a batch of knockout combos."""
    return [run_fba_knockout(model, combo) for combo in combos]


# ============================================================
# 2. XGBoost Surrogate
# ============================================================

def train_xgb(X, y, n_estimators=200, max_depth=6, lr=0.1):
    """Train XGBoost regressor."""
    m = xgb.XGBRegressor(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=lr,
        objective='reg:squarederror',
        random_state=42,
        verbosity=0,
    )
    m.fit(X, y)
    return m

def train_xgb_ensemble(X, y, n_models=5, n_estimators=100, max_depth=6, lr=0.1):
    """Train ensemble of XGBoost models with different seeds."""
    models = []
    for seed in range(n_models):
        m = xgb.XGBRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            learning_rate=lr,
            objective='reg:squarederror',
            random_state=seed * 100,
            verbosity=0,
        )
        m.fit(X, y)
        models.append(m)
    return models

def ensemble_predict(models, X):
    """Predict with ensemble, return (mean, std)."""
    preds = np.array([m.predict(X) for m in models])
    return preds.mean(axis=0), preds.std(axis=0)

def evaluate_r2(y_true, y_pred):
    """Compute R2 score."""
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    if ss_tot == 0:
        return 0.0
    return 1 - ss_res / ss_tot


# ============================================================
# 3. Input-space Active Learning Strategies
# ============================================================

def hamming_diversity_select(selected_masks, candidate_masks, n_select):
    """Select candidates with maximum Hamming distance from already selected."""
    if len(selected_masks) == 0:
        indices = np.random.choice(len(candidate_masks), n_select, replace=False)
        return indices

    dist_matrix = cdist(candidate_masks, selected_masks, metric='hamming')
    min_dists = dist_matrix.min(axis=1)  # distance to nearest selected

    top_indices = np.argsort(min_dists)[-n_select:]
    return top_indices

def ucb_select(models, candidate_masks, n_select, beta=2.0):
    """UCB-based selection: exploitation + exploration."""
    mean, std = ensemble_predict(models, candidate_masks)
    mean_norm = (mean - mean.min()) / (mean.max() - mean.min() + 1e-8)
    std_norm = (std - std.min()) / (std.max() - std.min() + 1e-8)
    ucb_score = mean_norm + beta * std_norm
    top_indices = np.argsort(ucb_score)[-n_select:]
    return top_indices


# ============================================================
# 4. Main Experiment
# ============================================================

def run_experiment():
    print("=" * 70)
    print("Step 1: Input-space Active Learning Micro-Validation")
    print("=" * 70)

    model = load_model()
    wt_growth = model.optimize().objective_value
    print(f"[Model] textbook loaded, WT growth = {wt_growth:.4f}")

    # --- Generate candidate pool ---
    N_POOL = 2000
    N_INITIAL = 200
    N_ROUNDS = 5
    N_SELECT = 20

    print(f"\n[Data] Generating {N_POOL} candidate knockout combos...")
    t0 = time.time()
    pool_combos = generate_random_combos(model, N_POOL, min_k=1, max_k=5, seed=42)
    pool_masks, gene_ids = build_knockout_mask(model, pool_combos)
    print(f"  Pool mask shape: {pool_masks.shape}, time: {time.time()-t0:.1f}s")

    # --- Compute FBA for all ---
    print(f"\n[FBA] Computing growth rates for {N_POOL} combos...")
    t0 = time.time()
    pool_growth = np.array(compute_growth_batch(model, pool_combos), dtype=np.float32)
    print(f"  Done in {time.time()-t0:.1f}s")
    print(f"  Growth stats: mean={pool_growth.mean():.4f}, std={pool_growth.std():.4f}")
    print(f"  Feasible: {(pool_growth > 1e-6).sum()}/{N_POOL}")

    # --- Run AL vs Random comparison ---
    results = {"al": [], "random": []}

    for strategy_name in ["random", "al"]:
        print(f"\n{'='*50}")
        print(f"Strategy: {strategy_name.upper()}")
        print(f"{'='*50}")

        rng = np.random.RandomState(42)
        initial_idx = rng.choice(N_POOL, N_INITIAL, replace=False)
        selected_idx = list(initial_idx)

        X_train = pool_masks[selected_idx]
        y_train = pool_growth[selected_idx]

        xgb_model = train_xgb(X_train, y_train)
        y_pred = xgb_model.predict(pool_masks)
        r2 = evaluate_r2(pool_growth, y_pred)
        print(f"  Round 0: {len(selected_idx)} samples, R2 = {r2:.4f}")

        round_results = [{"round": 0, "n_samples": len(selected_idx), "r2": r2}]

        for round_i in range(1, N_ROUNDS + 1):
            candidate_idx = [i for i in range(N_POOL) if i not in selected_idx]
            candidate_masks = pool_masks[candidate_idx]

            if strategy_name == "random":
                new_local_idx = rng.choice(len(candidate_idx), N_SELECT, replace=False)
            else:
                selected_masks = pool_masks[selected_idx]
                if r2 < 0.3:
                    new_local_idx = hamming_diversity_select(
                        selected_masks, candidate_masks, N_SELECT)
                else:
                    ensemble = train_xgb_ensemble(X_train, y_train, n_models=5)
                    new_local_idx = ucb_select(
                        ensemble, candidate_masks, N_SELECT, beta=2.0)

            new_global_idx = [candidate_idx[i] for i in new_local_idx]
            selected_idx.extend(new_global_idx)

            X_train = pool_masks[selected_idx]
            y_train = pool_growth[selected_idx]

            xgb_model = train_xgb(X_train, y_train)
            y_pred = xgb_model.predict(pool_masks)
            r2 = evaluate_r2(pool_growth, y_pred)

            print(f"  Round {round_i}: {len(selected_idx)} samples, R2 = {r2:.4f}")
            round_results.append({
                "round": round_i,
                "n_samples": len(selected_idx),
                "r2": float(r2),
            })

        results[strategy_name] = round_results

    # --- Summary ---
    print(f"\n{'='*70}")
    print("SUMMARY: Input-space AL vs Random")
    print(f"{'='*70}")
    print(f"{'Round':<8} {'N_samples':<12} {'AL R2':<12} {'Random R2':<12} {'Delta':<12}")
    print("-" * 56)
    for al_r, rand_r in zip(results["al"], results["random"]):
        delta = al_r["r2"] - rand_r["r2"]
        print(f"{al_r['round']:<8} {al_r['n_samples']:<12} "
              f"{al_r['r2']:<12.4f} {rand_r['r2']:<12.4f} {delta:<12.4f}")

    # --- Learning curve ---
    print(f"\n{'='*70}")
    print("LEARNING CURVE: XGBoost R2 vs Training Set Size")
    print(f"{'='*70}")
    sample_sizes = [50, 100, 200, 500, 1000, 2000]
    rng2 = np.random.RandomState(99)
    for n in sample_sizes:
        idx = rng2.choice(N_POOL, min(n, N_POOL), replace=False)
        X = pool_masks[idx]
        y = pool_growth[idx]
        xgb_m = train_xgb(X, y)
        y_pred = xgb_m.predict(pool_masks)
        r2 = evaluate_r2(pool_growth, y_pred)
        print(f"  N={n:>5}: R2 = {r2:.4f}")

    # --- Feature importance ---
    print(f"\n{'='*70}")
    print("FEATURE IMPORTANCE (XGBoost, full pool)")
    print(f"{'='*70}")
    xgb_full = train_xgb(pool_masks, pool_growth)
    importance = xgb_full.feature_importances_
    top_k = 20
    top_indices = np.argsort(importance)[-top_k:][::-1]
    print(f"  Top {top_k} genes by importance:")
    for rank, idx in enumerate(top_indices):
        print(f"    {rank+1:>3}. gene={gene_ids[idx]:<15} importance={importance[idx]:.4f}")

    # --- Hamming distance analysis ---
    print(f"\n{'='*70}")
    print("HAMMING DISTANCE ANALYSIS")
    print(f"{'='*70}")
    # Compare diversity of AL-selected vs random-selected
    rng3 = np.random.RandomState(42)
    al_idx_sample = rng3.choice(N_POOL, 300, replace=False)
    rand_idx_sample = rng3.choice(N_POOL, 300, replace=False)

    # Within-set pairwise Hamming distances (subsample for speed)
    n_sub = 50
    al_sub = pool_masks[al_idx_sample[:n_sub]]
    rand_sub = pool_masks[rand_idx_sample[:n_sub]]

    al_hamming = cdist(al_sub, al_sub, metric='hamming')
    rand_hamming = cdist(rand_sub, rand_sub, metric='hamming')

    # Exclude diagonal
    np.fill_diagonal(al_hamming, np.nan)
    np.fill_diagonal(rand_hamming, np.nan)

    print(f"  Mean pairwise Hamming (random subset): {np.nanmean(rand_hamming):.4f}")
    print(f"  Mean pairwise Hamming (AL-style subset): {np.nanmean(al_hamming):.4f}")

    return results


if __name__ == "__main__":
    results = run_experiment()
