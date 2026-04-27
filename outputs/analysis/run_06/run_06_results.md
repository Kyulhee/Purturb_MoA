# Run 06 Results: Norman 2019 Real Data Ablation Experiments

**Date:** 2026-04-27 | **Status:** Completed

---

## Summary

6 ablation configs tested on Norman 2019 (89K cells, 500 HVGs, 128 double-KO pairs).

| # | Config | best_corr | best_R2 | corr_add | corr_mul | R2_add | R2_mul |
|---|--------|-----------|---------|----------|----------|--------|--------|
| 1 | FCR baseline (no ICM) | 0.9552 | 0.8811 | 0.9547 | 0.9459 | 0.8806 | 0.8515 |
| 2 | FCR + ICM | 0.9538 | 0.8846 | 0.9528 | 0.9448 | 0.8827 | 0.8645 |
| 3 | FCR + linear z_tx head | **0.9571** | 0.8844 | **0.9566** | 0.9451 | 0.8827 | 0.8545 |
| 4 | FCR + ICM + linear z_tx | 0.9561 | **0.8861** | 0.9557 | 0.9396 | **0.8855** | 0.8508 |
| 5 | FCR + ICM + comp loss | 0.9464 | 0.8631 | 0.9456 | 0.9349 | 0.8602 | 0.8288 |
| 6 | ICM + linear + comp loss | 0.9550 | 0.8740 | 0.9545 | 0.9397 | 0.8729 | 0.8288 |

## Key Findings

### 1. All configs perform similarly on real data (R2 range: 0.86-0.89)

Unlike synthetic data where comp loss caused a dramatic improvement (0.20 -> 0.79), on real data all configs achieve R2 ~0.88. The baseline FCR is already near ceiling for gene-space composition.

### 2. Comp consistency loss HURTS performance on real data

Config 5 (ICM + comp loss) gets the worst performance (R2=0.8631 vs baseline 0.8811). The explicit composition constraint conflicts with the encoder-decoder's naturally learned representation.

**Why the synthetic-real discrepancy?**
- Synthetic: ground truth z_tx has known composition rules, but encoder breaks them → comp loss enforces the rules → helps
- Real: no ground truth z_tx; encoder-decoder jointly learn to compose in gene space → comp loss over-constrains the representation → hurts

### 3. Linear z_tx head provides marginal benefit

Config 3-4 (linear z_tx head) slightly outperform on correlation (0.957 vs 0.955). But the difference is within noise.

### 4. ICM provides marginal benefit on single cell type

With only 1 cell type, ICM acts as a structural regularizer (z_tx variance constraint). Slight R2 improvement (0.881 -> 0.885).

### 5. Additive composition dominates (110/128 pairs)

Additive composition (z_tx_1 + z_tx_2) is better than multiplicative for 110 out of 128 double-KO pairs (86%). This is consistent with the biological expectation that independent perturbations combine additively in expression space.

### Top performing pairs (config 4):
- MAP2K6+ELMSAN1: R2=0.945
- BCL2L11+BAK1: R2=0.943
- MAP2K3+MAP2K6: R2=0.944
- KLF1+MAP2K6: R2=0.941

## Implications

1. **Gene-space composition is robust**: The encoder-decoder naturally learns to compose. No explicit composition regularization needed on real data.
2. **Comp loss is only useful for synthetic diagnostics**: Where ground truth z_tx exists and we want to enforce latent-space compositionality.
3. **The RQ2 question is answered positively**: Composition works on real data (R2=0.88) without special mechanisms.
4. **For the paper**: FCR baseline is sufficient for compositionality. ICM is critical for invariance/transfer (RQ1/RQ3), not composition (RQ2).
