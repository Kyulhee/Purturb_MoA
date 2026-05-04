# BioEval Run 20d — GEARS DL Model Analysis

**Date**: 2026-05-03

## K562

**Perturbations**: 1092, **Genes**: 5000

| Model | R2 | Pearson | Dir_deg | DEG_auprc | f1@50 | MSE |
|-------|----|---------|---------|-----------|-------|-----|
| GEARS | 0.085 | 0.513 | 0.888 | 0.345 | 0.241 | 0.0057 |
| Ridge | 0.610 | 0.771 | 0.988 | 0.633 | 0.380 | 0.0022 |
| mean_predictor | -0.027 | 0.000 | 0.000 | 0.024 | 0.007 | 0.0082 |
| mean_effect | 0.038 | 0.384 | 0.789 | 0.186 | 0.150 | 0.0070 |

### H3: GEARS vs Baselines (K562)

- vs Ridge: GEARS wins 0/4, loses 4/4
- vs mean_predictor: GEARS wins 4/4, loses 0/4
- vs mean_effect: GEARS wins 4/4, loses 0/4

### H1: MSE vs Dir_deg Ranking (K562)

- MSE rank: {'GEARS': 2, 'Ridge': 1, 'mean_effect': 3, 'mean_predictor': 4}
- Dir_deg rank: {'GEARS': 2, 'Ridge': 1, 'mean_effect': 3, 'mean_predictor': 4}
- Kendall τ(MSE_rank, Dir_deg_rank) = 1.000, p = 0.083
- No significant rank reversal (|τ| >= 0.7)

## RPE1

**Perturbations**: 1543, **Genes**: 5000

| Model | R2 | Pearson | Dir_deg | DEG_auprc | f1@50 | MSE |
|-------|----|---------|---------|-----------|-------|-----|
| GEARS | 0.147 | 0.613 | 0.890 | 0.481 | 0.193 | 0.0100 |
| Ridge | 0.696 | 0.824 | 0.983 | 0.724 | 0.281 | 0.0036 |
| mean_predictor | -0.027 | 0.000 | 0.000 | 0.066 | 0.007 | 0.0201 |
| mean_effect | -0.148 | 0.627 | 0.874 | 0.517 | 0.205 | 0.0118 |

### H3: GEARS vs Baselines (RPE1)

- vs Ridge: GEARS wins 0/4, loses 4/4
- vs mean_predictor: GEARS wins 4/4, loses 0/4
- vs mean_effect: GEARS wins 2/4, loses 2/4

### H1: MSE vs Dir_deg Ranking (RPE1)

- MSE rank: {'GEARS': 2, 'Ridge': 1, 'mean_effect': 3, 'mean_predictor': 4}
- Dir_deg rank: {'GEARS': 2, 'Ridge': 1, 'mean_effect': 3, 'mean_predictor': 4}
- Kendall τ(MSE_rank, Dir_deg_rank) = 1.000, p = 0.083
- No significant rank reversal (|τ| >= 0.7)

## Norman

**Perturbations**: 283, **Genes**: 5045

| Model | R2 | Pearson | Dir_deg | DEG_auprc | f1@50 | MSE |
|-------|----|---------|---------|-----------|-------|-----|
| GEARS | -0.699 | 0.256 | 0.422 | 0.182 | 0.154 | 0.0074 |
| Ridge | 0.896 | 0.962 | 1.000 | 0.925 | 0.551 | 0.0004 |
| mean_predictor | -0.004 | 0.000 | 0.000 | 0.017 | 0.009 | 0.0053 |
| mean_effect | -0.123 | 0.509 | 0.833 | 0.423 | 0.334 | 0.0040 |

### H3: GEARS vs Baselines (Norman)

- vs Ridge: GEARS wins 0/4, loses 4/4
- vs mean_predictor: GEARS wins 3/4, loses 1/4
- vs mean_effect: GEARS wins 0/4, loses 4/4

### H1: MSE vs Dir_deg Ranking (Norman)

- MSE rank: {'GEARS': 4, 'Ridge': 1, 'mean_effect': 2, 'mean_predictor': 3}
- Dir_deg rank: {'GEARS': 3, 'Ridge': 1, 'mean_effect': 2, 'mean_predictor': 4}
- Kendall τ(MSE_rank, Dir_deg_rank) = 0.667, p = 0.333
- **RANK REVERSAL DETECTED** (|τ| < 0.7)

---

## Overall H3 Summary

- GEARS wins 17/36 comparisons across all datasets and baselines
- GEARS datasets available: K562, RPE1, Norman
- Ridge wins 12/12 vs mean_predictor (confirmed from run_16-19)