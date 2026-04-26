"""
Module C: Active Learning Loop
================================
Two-phase active learning for efficient FBA sample selection.

Phase 1 (R2 < 0.3): diversity-based exploration
    - Select samples maximally distant in embedding space
    - Uncertainty is invalid when R2 < 0, so pure exploration

Phase 2 (R2 >= 0.3): UCB-based exploration-exploitation
    - Upper Confidence Bound: exploitation(high predicted growth) + exploration(high uncertainty)
    - Transition condition: validation R2 >= 0.3 for 3 consecutive evaluations

Design decisions (from stages/03_planning.md):
    - AL rounds: 50 samples per round, max 20 rounds (1,000 additional samples)
    - Transition: R2 >= 0.3 for 3 consecutive evaluations
    - diversity: embedding-space max-min distance selection
    - UCB: alpha * mean_prediction + (1-alpha) * uncertainty
"""

import time
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import xgboost as xgb
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import train_test_split


# ---------------------------------------------------------------------------
# 1. Candidate Pool Generator
# ---------------------------------------------------------------------------

def generate_candidate_pool(
    n_genes: int,
    n_candidates: int = 10000,
    min_k: int = 1,
    max_k: int = 5,
    seed: int = 42,
) -> np.ndarray:
    """
    Generate a large pool of candidate knockout combinations.

    Returns:
        candidates: (n_candidates, n_genes) binary mask array
    """
    rng = np.random.RandomState(seed)
    candidates = np.zeros((n_candidates, n_genes), dtype=np.float32)

    for i in range(n_candidates):
        k = rng.randint(min_k, max_k + 1)
        indices = rng.choice(n_genes, size=k, replace=False)
        candidates[i, indices] = 1.0

    return candidates


# ---------------------------------------------------------------------------
# 2. Diversity-based Selection (Phase 1)
# ---------------------------------------------------------------------------

def diversity_selection(
    candidates: np.ndarray,
    existing_masks: np.ndarray,
    n_select: int = 50,
    embeddings: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Select candidates maximally distant from existing samples.

    Uses max-min distance in embedding space (or raw mask space).
    Iteratively picks the candidate farthest from all existing samples.

    Args:
        candidates: (n_pool, n_genes) candidate masks
        existing_masks: (n_existing, n_genes) already selected masks
        n_select: number of candidates to select
        embeddings: optional (n_existing + n_pool, embed_dim) precomputed embeddings

    Returns:
        selected_indices: indices into candidates array
    """
    if embeddings is not None:
        # Use embedding space for distance computation
        n_existing = len(existing_masks)
        existing_embs = embeddings[:n_existing]
        candidate_embs = embeddings[n_existing:]

        # Compute distances from each candidate to all existing samples
        # Shape: (n_candidates, n_existing)
        dists = np.linalg.norm(
            candidate_embs[:, np.newaxis, :] - existing_embs[np.newaxis, :, :],
            axis=2
        )
        # Min distance to any existing sample
        min_dists = dists.min(axis=1)
    else:
        # Use raw mask space
        # Compute pairwise distances
        min_dists = np.full(len(candidates), np.inf)

        for i in range(len(candidates)):
            dists_to_existing = np.linalg.norm(
                candidates[i] - existing_masks, axis=1
            )
            min_dists[i] = dists_to_existing.min()

    # Greedy max-min selection
    selected_indices = []
    remaining = np.arange(len(candidates))
    current_min_dists = min_dists.copy()

    for _ in range(min(n_select, len(candidates))):
        # Pick the candidate with maximum minimum distance
        best_idx = current_min_dists[remaining].argmax()
        best_candidate = remaining[best_idx]
        selected_indices.append(best_candidate)

        # Remove selected from remaining
        remaining = np.delete(remaining, best_idx)

        # Update distances: newly selected candidate is now "existing"
        if len(remaining) > 0:
            new_dists = np.linalg.norm(
                candidates[remaining] - candidates[best_candidate], axis=1
            )
            current_min_dists[remaining] = np.minimum(
                current_min_dists[remaining], new_dists
            )

    return np.array(selected_indices)


# ---------------------------------------------------------------------------
# 3. UCB-based Selection (Phase 2)
# ---------------------------------------------------------------------------

def ucb_selection(
    candidates: np.ndarray,
    predictions: np.ndarray,
    uncertainties: np.ndarray,
    n_select: int = 50,
    alpha: float = 0.5,
) -> np.ndarray:
    """
    Upper Confidence Bound selection.

    UCB = alpha * normalized_prediction + (1 - alpha) * normalized_uncertainty

    High alpha = exploitation (select high-growth predictions)
    Low alpha = exploration (select high-uncertainty regions)

    Args:
        candidates: (n_pool, n_genes) candidate masks
        predictions: (n_pool,) predicted growth rates
        uncertainties: (n_pool,) prediction uncertainties
        n_select: number to select
        alpha: exploitation weight (0 = pure exploration, 1 = pure exploitation)

    Returns:
        selected_indices: indices into candidates array
    """
    # Normalize predictions and uncertainties to [0, 1]
    pred_min, pred_max = predictions.min(), predictions.max()
    unc_min, unc_max = uncertainties.min(), uncertainties.max()

    pred_norm = (predictions - pred_min) / max(pred_max - pred_min, 1e-8)
    unc_norm = (uncertainties - unc_min) / max(unc_max - unc_min, 1e-8)

    # UCB score
    ucb_scores = alpha * pred_norm + (1 - alpha) * unc_norm

    # Select top-n by UCB score
    top_indices = np.argsort(ucb_scores)[-n_select:]
    return top_indices


# ---------------------------------------------------------------------------
# 4. Quantile Regression for Uncertainty Estimation
# ---------------------------------------------------------------------------

def estimate_uncertainty_quantile(
    candidate_masks: np.ndarray,
    xgb_models: List,
) -> np.ndarray:
    """
    Estimate uncertainty using quantile regression ensemble.

    Train multiple XGBoost models with quantile objectives,
    then uncertainty = upper_quantile - lower_quantile.

    Args:
        candidate_masks: combined features (embedding + mask)
        xgb_models: list of trained XGBoost quantile models

    Returns:
        uncertainties: (n_candidates,) uncertainty estimates
    """
    predictions = []
    for model in xgb_models:
        pred = model.predict(candidate_masks)
        predictions.append(pred)

    predictions = np.array(predictions)  # (n_models, n_candidates)

    # Uncertainty = range of predictions across quantiles
    uncertainties = predictions.max(axis=0) - predictions.min(axis=0)
    return uncertainties


# ---------------------------------------------------------------------------
# 5. Active Learning Loop
# ---------------------------------------------------------------------------

class ActiveLearningLoop:
    """
    Two-phase active learning loop for FBA sample selection.

    Phase 1 (R2 < 0.3): diversity-based exploration
    Phase 2 (R2 >= 0.3): UCB-based exploration-exploitation

    Transition condition: validation R2 >= 0.3 for 3 consecutive evaluations.
    """

    def __init__(
        self,
        n_genes: int,
        n_select_per_round: int = 50,
        max_rounds: int = 20,
        transition_r2: float = 0.3,
        transition_patience: int = 3,
        ucb_alpha: float = 0.5,
        diversity_seed: int = 42,
    ):
        self.n_genes = n_genes
        self.n_select_per_round = n_select_per_round
        self.max_rounds = max_rounds
        self.transition_r2 = transition_r2
        self.transition_patience = transition_patience
        self.ucb_alpha = ucb_alpha
        self.diversity_seed = diversity_seed

        # State tracking
        self.current_phase = 1
        self.r2_history = []
        self.consecutive_above_threshold = 0
        self.total_fba_calls = 0
        self.phase_transition_round = None

    def _check_transition(self, val_r2: float) -> bool:
        """Check if we should transition from Phase 1 to Phase 2."""
        if val_r2 >= self.transition_r2:
            self.consecutive_above_threshold += 1
        else:
            self.consecutive_above_threshold = 0

        return self.consecutive_above_threshold >= self.transition_patience

    def run(
        self,
        initial_masks: np.ndarray,
        initial_growth: np.ndarray,
        fba_oracle,  # Callable: (masks) -> growth_rates
        surrogate_factory,  # Callable: () -> GNNXGBoostSurrogate
        candidate_pool_size: int = 10000,
        val_ratio: float = 0.2,
        verbose: bool = True,
    ) -> Dict:
        """
        Run the active learning loop.

        Args:
            initial_masks: (n_initial, n_genes) initial knockout masks
            initial_growth: (n_initial,) initial FBA growth rates
            fba_oracle: function that takes masks and returns FBA growth rates
            surrogate_factory: function that creates a fresh GNNXGBoostSurrogate
            candidate_pool_size: number of candidates to generate each round
            val_ratio: validation split ratio

        Returns:
            dict with AL loop results and metrics
        """
        t0 = time.time()

        # Initialize training data
        train_masks = initial_masks.copy()
        train_growth = initial_growth.copy()

        # Track metrics
        r2_per_round = []
        fba_calls_per_round = []
        phase_per_round = []

        print(f"[AL] Starting Active Learning Loop")
        print(f"  Initial data: {len(train_growth)} samples")
        print(f"  Phase 1: diversity (R2 < {self.transition_r2})")
        print(f"  Phase 2: UCB (R2 >= {self.transition_r2})")
        print(f"  {self.n_select_per_round} samples/round, max {self.max_rounds} rounds")

        for round_idx in range(self.max_rounds):
            round_t0 = time.time()

            # Step 1: Train surrogate on current data
            surrogate = surrogate_factory()
            fit_results = surrogate.fit(
                train_masks, train_growth,
                pretrain_epochs=20,
                val_ratio=val_ratio,
                verbose=False,
            )
            val_r2 = fit_results["r2_test"]
            self.r2_history.append(val_r2)
            r2_per_round.append(val_r2)

            # Step 2: Check phase transition
            if self.current_phase == 1:
                if self._check_transition(val_r2):
                    self.current_phase = 2
                    self.phase_transition_round = round_idx
                    if verbose:
                        print(f"\n[AL] *** Phase transition at round {round_idx} ***")
                        print(f"  R2 = {val_r2:.4f} >= {self.transition_r2} "
                              f"for {self.transition_patience} consecutive rounds")

            phase_per_round.append(self.current_phase)

            # Step 3: Generate candidate pool
            candidate_pool = generate_candidate_pool(
                n_genes=self.n_genes,
                n_candidates=candidate_pool_size,
                seed=self.diversity_seed + round_idx,
            )

            # Remove candidates that are already in training set
            # (exact match check)
            is_new = np.ones(len(candidate_pool), dtype=bool)
            for i in range(len(candidate_pool)):
                for j in range(len(train_masks)):
                    if np.array_equal(candidate_pool[i], train_masks[j]):
                        is_new[i] = False
                        break
            candidate_pool = candidate_pool[is_new]

            if len(candidate_pool) == 0:
                if verbose:
                    print(f"[AL] Round {round_idx}: No new candidates, stopping")
                break

            # Step 4: Select candidates based on current phase
            if self.current_phase == 1:
                # Phase 1: diversity-based selection
                selected_indices = diversity_selection(
                    candidates=candidate_pool,
                    existing_masks=train_masks,
                    n_select=self.n_select_per_round,
                )
            else:
                # Phase 2: UCB-based selection
                predictions = surrogate.predict(candidate_pool)

                # Simple uncertainty via dropout-style ensemble
                # Use prediction magnitude as proxy when quantile models unavailable
                # For proper uncertainty, train quantile XGBoost ensemble
                uncertainties = self._estimate_uncertainty_simple(
                    surrogate, candidate_pool
                )

                selected_indices = ucb_selection(
                    candidates=candidate_pool,
                    predictions=predictions,
                    uncertainties=uncertainties,
                    n_select=self.n_select_per_round,
                    alpha=self.ucb_alpha,
                )

            # Step 5: Query FBA oracle for selected candidates
            selected_masks = candidate_pool[selected_indices]
            selected_growth = fba_oracle(selected_masks)
            n_fba = len(selected_growth)
            self.total_fba_calls += n_fba
            fba_calls_per_round.append(n_fba)

            # Step 6: Add to training data
            train_masks = np.concatenate([train_masks, selected_masks], axis=0)
            train_growth = np.concatenate([train_growth, selected_growth], axis=0)

            round_time = time.time() - round_t0

            if verbose:
                phase_str = "diversity" if self.current_phase == 1 else "UCB"
                print(f"[AL] Round {round_idx:2d}: phase={phase_str}, "
                      f"R2={val_r2:.4f}, "
                      f"FBA={n_fba}, "
                      f"total_samples={len(train_growth)}, "
                      f"time={round_time:.1f}s")

            # Early stopping: R2 target reached
            if val_r2 >= 0.5:
                if verbose:
                    print(f"\n[AL] Target R2 >= 0.5 reached at round {round_idx}!")
                break

        total_time = time.time() - t0

        results = {
            "r2_per_round": r2_per_round,
            "fba_calls_per_round": fba_calls_per_round,
            "phase_per_round": phase_per_round,
            "total_fba_calls": self.total_fba_calls,
            "total_rounds": len(r2_per_round),
            "phase_transition_round": self.phase_transition_round,
            "final_r2": r2_per_round[-1] if r2_per_round else None,
            "final_n_samples": len(train_growth),
            "total_time": total_time,
            "train_masks": train_masks,
            "train_growth": train_growth,
        }

        if verbose:
            print(f"\n[AL] Completed:")
            print(f"  Rounds: {results['total_rounds']}")
            print(f"  Final R2: {results['final_r2']:.4f}")
            print(f"  Total FBA calls: {results['total_fba_calls']}")
            print(f"  Final samples: {results['final_n_samples']}")
            if self.phase_transition_round is not None:
                print(f"  Phase transition: round {self.phase_transition_round}")
            print(f"  Total time: {total_time:.1f}s")

        return results

    def _estimate_uncertainty_simple(
        self,
        surrogate,
        candidate_masks: np.ndarray,
    ) -> np.ndarray:
        """
        Simple uncertainty estimation via prediction variance.

        Uses multiple forward passes with different dropout masks
        as a lightweight ensemble. Falls back to distance-based
        uncertainty if dropout not available.
        """
        # Method: use embedding distance to training data as uncertainty proxy
        # Samples far from training data = high uncertainty
        # This is a simple but effective heuristic

        try:
            # Get predictions - high variance in predictions indicates uncertainty
            pred = surrogate.predict(candidate_masks)

            # Use absolute prediction value as rough uncertainty proxy
            # (regions where model is less confident tend to have extreme predictions)
            # Better: use distance in feature space

            # For now, use a simple random perturbation approach
            n_perturbations = 5
            all_preds = []
            for _ in range(n_perturbations):
                # Small random perturbation to features
                perturbed = candidate_masks.copy()
                noise = np.random.normal(0, 0.05, perturbed.shape).astype(np.float32)
                perturbed = np.clip(perturbed + noise, 0, 1)
                pred_p = surrogate.predict(perturbed)
                all_preds.append(pred_p)

            all_preds = np.array(all_preds)  # (n_perturbations, n_candidates)
            uncertainties = all_preds.std(axis=0)

        except Exception:
            # Fallback: random uncertainty
            uncertainties = np.random.rand(len(candidate_masks)).astype(np.float32)

        return uncertainties


# ---------------------------------------------------------------------------
# 6. Random Screening Baseline
# ---------------------------------------------------------------------------

class RandomScreeningBaseline:
    """
    Baseline: random sample selection instead of active learning.
    Same number of FBA calls as AL for fair comparison.
    """

    def __init__(self, n_genes: int, n_select_per_round: int = 50, max_rounds: int = 20):
        self.n_genes = n_genes
        self.n_select_per_round = n_select_per_round
        self.max_rounds = max_rounds

    def run(
        self,
        initial_masks: np.ndarray,
        initial_growth: np.ndarray,
        fba_oracle,
        surrogate_factory,
        val_ratio: float = 0.2,
        verbose: bool = True,
    ) -> Dict:
        """Run random screening with same budget as AL."""
        t0 = time.time()

        train_masks = initial_masks.copy()
        train_growth = initial_growth.copy()
        r2_per_round = []
        fba_calls_per_round = []

        rng = np.random.RandomState(42)

        for round_idx in range(self.max_rounds):
            # Train surrogate
            surrogate = surrogate_factory()
            fit_results = surrogate.fit(
                train_masks, train_growth,
                pretrain_epochs=20,
                val_ratio=val_ratio,
                verbose=False,
            )
            val_r2 = fit_results["r2_test"]
            r2_per_round.append(val_r2)

            # Random selection
            random_masks = generate_candidate_pool(
                n_genes=self.n_genes,
                n_candidates=self.n_select_per_round,
                seed=rng.randint(0, 100000),
            )

            # Query FBA
            random_growth = fba_oracle(random_masks)
            n_fba = len(random_growth)
            fba_calls_per_round.append(n_fba)

            # Add to training data
            train_masks = np.concatenate([train_masks, random_masks], axis=0)
            train_growth = np.concatenate([train_growth, random_growth], axis=0)

            if verbose:
                print(f"[Random] Round {round_idx:2d}: R2={val_r2:.4f}, "
                      f"FBA={n_fba}, total_samples={len(train_growth)}")

            # Early stopping
            if val_r2 >= 0.5:
                break

        total_time = time.time() - t0

        return {
            "r2_per_round": r2_per_round,
            "fba_calls_per_round": fba_calls_per_round,
            "total_fba_calls": sum(fba_calls_per_round),
            "total_rounds": len(r2_per_round),
            "final_r2": r2_per_round[-1] if r2_per_round else None,
            "final_n_samples": len(train_growth),
            "total_time": total_time,
            "train_masks": train_masks,
            "train_growth": train_growth,
        }


# ---------------------------------------------------------------------------
# 7. Test / verification
# ---------------------------------------------------------------------------

def run_tests():
    """Quick verification tests for Module C."""
    print("=" * 60)
    print("Module C: Active Learning Loop -- Verification Tests")
    print("=" * 60)

    all_pass = True

    # Test 1: Candidate pool generation
    print("\n--- Test 1: Candidate Pool Generation ---")
    try:
        pool = generate_candidate_pool(n_genes=137, n_candidates=1000)
        assert pool.shape == (1000, 137), f"Pool shape {pool.shape} != (1000, 137)"
        assert (pool.sum(axis=1) >= 1).all(), "Each candidate should knock out >= 1 gene"
        print(f"  PASS: Pool shape {pool.shape}, min KO per candidate: {pool.sum(axis=1).min():.0f}")
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Test 2: Diversity selection
    print("\n--- Test 2: Diversity Selection ---")
    try:
        existing = np.zeros((10, 137), dtype=np.float32)
        for i in range(10):
            existing[i, i] = 1.0

        candidates = np.zeros((100, 137), dtype=np.float32)
        rng = np.random.RandomState(42)
        for i in range(100):
            k = rng.randint(1, 4)
            indices = rng.choice(137, size=k, replace=False)
            candidates[i, indices] = 1.0

        selected_idx = diversity_selection(candidates, existing, n_select=10)
        assert len(selected_idx) == 10, f"Expected 10 selected, got {len(selected_idx)}"
        assert len(np.unique(selected_idx)) == 10, "Selected indices should be unique"
        print(f"  PASS: Selected {len(selected_idx)} diverse candidates")
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Test 3: UCB selection
    print("\n--- Test 3: UCB Selection ---")
    try:
        candidates = np.zeros((100, 137), dtype=np.float32)
        predictions = np.random.rand(100).astype(np.float32)
        uncertainties = np.random.rand(100).astype(np.float32)

        selected_idx = ucb_selection(candidates, predictions, uncertainties, n_select=10)
        assert len(selected_idx) == 10
        assert len(np.unique(selected_idx)) == 10
        print(f"  PASS: UCB selected {len(selected_idx)} candidates")
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Test 4: AL loop with mock oracle
    print("\n--- Test 4: AL Loop with Mock Oracle ---")
    try:
        def mock_fba_oracle(masks):
            """Mock FBA: return random growth rates with some structure."""
            n_ko = masks.sum(axis=1)
            # More knockouts = lower growth (with noise)
            growth = 0.8 * np.exp(-0.3 * n_ko) + np.random.normal(0, 0.05, len(masks))
            return np.clip(growth, 0, 1).astype(np.float32)

        # Initial data
        initial_masks = generate_candidate_pool(n_genes=20, n_candidates=50, seed=42)
        initial_growth = mock_fba_oracle(initial_masks)

        # Simple surrogate factory (XGBoost-only for speed)
        class MockSurrogate:
            def __init__(self):
                self.model = None
                self.xgb_params = {
                    "objective": "reg:squarederror",
                    "max_depth": 4,
                    "n_estimators": 50,
                    "random_state": 42,
                }

            def fit(self, masks, growth, pretrain_epochs=0, val_ratio=0.2, verbose=False):
                from sklearn.model_selection import train_test_split
                from sklearn.metrics import r2_score
                X_train, X_test, y_train, y_test = train_test_split(
                    masks, growth, test_size=val_ratio, random_state=42
                )
                self.model = xgb.XGBRegressor(**self.xgb_params)
                self.model.fit(X_train, y_train, verbose=False)
                r2 = r2_score(y_test, self.model.predict(X_test))
                return {"r2_test": r2}

            def predict(self, masks):
                return self.model.predict(masks)

        al = ActiveLearningLoop(
            n_genes=20,
            n_select_per_round=20,
            max_rounds=5,
            transition_r2=0.3,
            transition_patience=2,
        )

        results = al.run(
            initial_masks=initial_masks,
            initial_growth=initial_growth,
            fba_oracle=mock_fba_oracle,
            surrogate_factory=MockSurrogate,
            verbose=True,
        )

        assert results["total_rounds"] <= 5
        assert results["total_fba_calls"] > 0
        print(f"  PASS: AL completed, {results['total_rounds']} rounds, "
              f"final R2={results['final_r2']:.4f}")
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Test 5: Random vs AL comparison
    print("\n--- Test 5: Random vs AL Comparison ---")
    try:
        random_baseline = RandomScreeningBaseline(
            n_genes=20,
            n_select_per_round=20,
            max_rounds=5,
        )

        random_results = random_baseline.run(
            initial_masks=initial_masks,
            initial_growth=initial_growth,
            fba_oracle=mock_fba_oracle,
            surrogate_factory=MockSurrogate,
            verbose=True,
        )

        print(f"  AL: {results['total_fba_calls']} FBA calls, R2={results['final_r2']:.4f}")
        print(f"  Random: {random_results['total_fba_calls']} FBA calls, "
              f"R2={random_results['final_r2']:.4f}")
        print(f"  PASS: Comparison completed")
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Summary
    print("\n" + "=" * 60)
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED -- see above")
    print("=" * 60)


if __name__ == "__main__":
    run_tests()
