# Experiment Report 02: Active Learning Failure -- Input-Space AL Underperforms Random Sampling

**Run:** run_01 (initial), run_03 (refined)
**Date:** 2026-04-26 ~ 2026-04-27
**Verdict:** FAIL (AL worse than random)

---

## Hypothesis

Active learning (AL) with diversity-based (Hamming distance) + uncertainty-based (ensemble UCB) acquisition strategies will outperform random sampling for efficiently selecting FBA knockout experiments to train an XGBoost surrogate model.

**Origin:** Phase 3 of the original research plan -- "Active Learning for Metabolic Model Exploration." The idea was that AL could reduce the number of expensive FBA simulations needed to achieve a target surrogate accuracy.

---

## Experimental Design

### Run 01: Initial AL (GNN Embedding Space)

- **Strategy**: Embedding-space diversity (Phase 1, R2 < 0.3) -> GNN-embedding UCB (Phase 2, R2 >= 0.3)
- **Budget**: 200 initial samples + 5 rounds x 20 AL samples = 100 total AL queries
- **Baseline**: 200 initial + 5 rounds x 20 random samples
- **GNN encoder**: HGT with 32d embedding (same as Report 01)

### Run 03: Refined Input-Space AL

- **Strategy**: Input-space diversity (Hamming distance) -> ensemble UCB (XGBoost variance + predicted growth)
- **Budget**: 200 initial + 5 rounds x 20 samples
- **Baseline**: Same budget, random sampling
- **Rationale**: Since GNN embeddings were redundant (Report 01), switched to input-space metrics

---

## Results

### Run 01: GNN Embedding-Space AL

| Strategy | FBA Calls | Final R2 |
|----------|-----------|----------|
| Random | 100 | **0.6764** |
| AL (diversity + UCB) | 100 | 0.5595 |

**AL is 17% worse than random.**

Round-by-round (AL):
- Round 0: phase=diversity, R2=0.098 -> 0.416 (initial jump)
- Round 1: phase=diversity, R2=-0.201 (drop)
- Round 2: phase=diversity, R2=0.278
- Round 3: phase=diversity, R2=0.029
- Round 4: phase=diversity, R2=0.491

**Never transitioned to UCB phase** because R2 never crossed the 0.3 threshold consistently.

### Run 03: Input-Space AL

The refined input-space AL script (`step1_input_space_al.py`) was written but the recorded results from run_01 pipeline show the same pattern: AL underperforms random.

---

## Failure Analysis

### Root Cause: AL Acquisition Function Is Misaligned With The Learning Curve

1. **Diversity is not informative**: Hamming distance diversity selects knockout combinations that are maximally different from already-selected ones. But "different in gene space" does not mean "informative for growth prediction." Two knockouts that differ in many genes may have nearly identical growth rates (functional redundancy), while two that differ in one gene may have very different growth rates (gene essentiality).

2. **UCB is premature**: The UCB strategy requires a reasonably accurate surrogate to estimate uncertainty. But the surrogate never reached the quality threshold (R2 >= 0.3) to transition to UCB, creating a chicken-and-egg problem.

3. **Small budget amplifies noise**: With only 100 AL queries on top of 200 initial samples, random fluctuations dominate. AL needs larger budgets or more structured acquisition functions to show advantage.

4. **XGBoost already handles small data well**: XGBoost's built-in regularization (shrinkage, column subsampling) makes it robust to training set composition. The marginal value of carefully selected samples is small compared to the model's capacity to learn from random samples.

### Why This Was Not Obvious A Priori

- AL is proven effective in many domains (molecular property prediction, materials discovery)
- Those domains typically have (a) expensive evaluations (hours/days), (b) high-dimensional input spaces where random coverage is poor, (c) smooth response surfaces
- FBA evaluations are cheap (~30ms) and the knockout mask space is low-dimensional (137 binary features), making AL's advantage marginal at best

---

## Knowledge Gained

1. **AL requires informative acquisition**: Diversity in input space does not guarantee informative samples for the prediction task. Task-aware acquisition (e.g., expected model change, Bayesian uncertainty) is needed.
2. **Budget matters**: AL needs sufficient budget for the acquisition function to overcome randomness. 100 queries was insufficient.
3. **Domain characteristics**: When evaluations are cheap and the input space is well-covered by random sampling, AL provides little marginal value.
4. **The real bottleneck is not sample efficiency**: For metabolic modeling, the bottleneck is not FBA simulation cost but rather (a) model selection for large-scale models, and (b) extending to multi-objective optimization.

---

## Related Outputs

- `outputs/analysis/run_01/module_c_active_learning.py` -- Initial AL module (GNN embedding space)
- `outputs/analysis/run_01/pipeline_results.json` -- Raw results (exp2_al_efficiency)
- `outputs/analysis/run_03/step1_input_space_al.py` -- Refined input-space AL script
