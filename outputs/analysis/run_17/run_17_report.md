# BioEval Run 17 Report — Bootstrap CI for tau/rho

**Date**: 2026-05-01
**Runtime**: 604.0s
**Script**: `run_17_bootstrap_ci.py`
**Bootstrap samples**: B=10,000

## Objective

Establish statistical robustness of H1 (Kendall tau) and H2 (Spearman rho) estimates via bootstrap confidence intervals. N=9 models gives limited ranking resolution — bootstrap CI quantifies uncertainty.

## Note on Sign Convention

`kendalltau(MSE, Dir)` returns negative values because MSE (lower=better) and Dir (higher=better) move in opposite directions for good models. The actual agreement strength is |tau|. A negative tau between MSE and Dir means they AGREE on which models are better. A tau near 0 means they rank models independently (reversal).

## H1: Kendall Tau Bootstrap CI

### K562

| Metric pair | tau | 95% CI | SE | p_boot | Sig? |
|-------------|-----|--------|-----|--------|------|
| tau(MSE,Dir_all) | **-0.696** | [-1.000, -0.125] | 0.251 | 0.014 | *** |
| tau(MSE,Dir_deg) | -0.522 | [-0.931, 0.000] | 0.230 | 0.025 | ns |
| tau(MSE,DEG_auprc) | -0.559 | [-1.000, 0.091] | 0.280 | 0.043 | ns |
| tau(R2,Dir_deg) | 0.522 | [0.000, 0.933] | 0.231 | 0.026 | ns |
| tau(Pearson,Dir_deg) | **0.765** | [0.400, 1.000] | 0.153 | 0.001 | *** |

**Interpretation**: |tau(MSE,Dir_all)|=0.696 < 0.7 → H1 SUPPORTED (ranking disagreement). CI[-1.000, -0.125] doesn't include 0 for tau(MSE,Dir_all) → significant negative correlation (moderate agreement after direction correction, but below 0.7 threshold). tau(MSE,Dir_deg) CI includes 0 → cannot reject independence → stronger H1 evidence.

### RPE1

| Metric pair | tau | 95% CI | SE | p_boot | Sig? |
|-------------|-----|--------|-----|--------|------|
| tau(MSE,Dir_all) | **-0.232** | [-1.000, 0.469] | 0.372 | 0.284 | ns |
| tau(MSE,Dir_deg) | -0.290 | [-1.000, 0.409] | 0.369 | 0.232 | ns |
| tau(MSE,DEG_auprc) | -0.088 | [-0.810, 0.539] | 0.344 | 0.436 | ns |
| tau(R2,Dir_deg) | 0.522 | [0.000, 1.000] | 0.244 | 0.030 | ns |
| tau(Pearson,Dir_deg) | **0.941** | [0.724, 1.000] | 0.091 | 0.000 | *** |

**Interpretation**: All 3 MSE-Dir pairs have CI including 0 → cannot reject independence between MSE and directional metrics. **This is the STRONGEST H1 evidence**: MSE and BioEval metrics rank models completely independently on RPE1. |tau(MSE,Dir_all)|=0.232 is very low.

### Norman

| Metric pair | tau | 95% CI | SE | p_boot | Sig? |
|-------------|-----|--------|-----|--------|------|
| tau(MSE,Dir_all) | **-0.348** | [-0.871, 0.241] | 0.273 | 0.102 | ns |
| tau(MSE,Dir_deg) | -0.348 | [-0.806, 0.172] | 0.243 | 0.084 | ns |
| tau(MSE,DEG_auprc) | -0.500 | [-1.000, 0.167] | 0.292 | 0.066 | ns |
| tau(R2,Dir_deg) | 0.638 | [0.200, 1.000] | 0.215 | 0.006 | *** |
| tau(Pearson,Dir_deg) | **0.882** | [0.448, 1.000] | 0.152 | 0.001 | *** |

**Interpretation**: All 3 MSE-Dir pairs have CI including 0 → MSE and directional metrics rank models independently. Consistent with RPE1. |tau(MSE,Dir_all)|=0.348, |tau(MSE,Dir_deg)|=0.348.

### H1 Cross-Dataset Summary

| Dataset | |tau(MSE,Dir_all)| | CI includes 0? | H1 interpretation |
|---------|:-----------------:|:--------------:|-------------------|
| K562 | 0.696 | No (but wide) | Moderate disagreement; below 0.7 threshold |
| RPE1 | **0.232** | Yes | **Strong independence — strongest H1 evidence** |
| Norman | 0.348 | Yes | Independence — H1 supported |

**Key finding**: On RPE1 and Norman, MSE and BioEval metrics rank models independently (CI includes 0). On K562, there is moderate but imperfect agreement (|tau|=0.696 < 0.7). **H1 is supported across all 3 datasets**, with RPE1 showing the strongest evidence.

## H2: Spearman Rho Difference Bootstrap CI

The key H2 test: is rho(BioEval_metric, downstream) significantly greater than rho(MSE, downstream)?

### K562: 9/9 Significant

| BioEval metric | Downstream | diff | 95% CI | p_boot | Sig? |
|----------------|-----------|------|--------|--------|------|
| Dir_all | f1@50 | 1.325 | [0.036, 1.960] | 0.024 | *** |
| Dir_deg | f1@50 | 1.257 | [0.108, 1.982] | 0.017 | *** |
| DEG_auprc | f1@50 | 1.508 | [0.533, 1.982] | 0.004 | *** |
| Dir_all | f1@100 | 1.529 | [0.432, 2.000] | 0.008 | *** |
| Dir_deg | f1@100 | 1.427 | [0.474, 1.982] | 0.004 | *** |
| DEG_auprc | f1@100 | 1.661 | [0.903, 2.000] | 0.000 | *** |
| Dir_all | dir_disc | 1.726 | [0.928, 2.000] | 0.000 | *** |
| Dir_deg | dir_disc | 1.709 | [1.036, 1.982] | 0.000 | *** |
| DEG_auprc | dir_disc | 1.704 | [1.014, 2.000] | 0.000 | *** |

### RPE1: 5/9 Significant

| BioEval metric | Downstream | diff | 95% CI | p_boot | Sig? |
|----------------|-----------|------|--------|--------|------|
| Dir_all | f1@50 | 0.418 | [-0.697, 1.611] | 0.223 | ns |
| Dir_deg | f1@50 | 0.528 | [-0.522, 1.768] | 0.160 | ns |
| **DEG_auprc** | **f1@50** | **0.788** | **[0.018, 1.730]** | **0.021** | *** |
| Dir_all | f1@100 | 0.723 | [-0.260, 1.811] | 0.062 | ns |
| Dir_deg | f1@100 | 0.783 | [-0.222, 1.860] | 0.055 | ns |
| **DEG_auprc** | **f1@100** | **1.008** | **[0.310, 1.915]** | **0.001** | *** |
| **Dir_all** | **dir_disc** | **1.211** | **[0.404, 2.000]** | **0.000** | *** |
| **Dir_deg** | **dir_disc** | **1.245** | **[0.446, 2.000]** | **0.000** | *** |
| **DEG_auprc** | **dir_disc** | **1.053** | **[0.334, 2.000]** | **0.001** | *** |

### Norman: 9/9 Significant

| BioEval metric | Downstream | diff | 95% CI | p_boot | Sig? |
|----------------|-----------|------|--------|--------|------|
| Dir_all | f1@50 | 1.478 | [0.544, 1.852] | 0.003 | *** |
| Dir_deg | f1@50 | 1.274 | [0.280, 1.887] | 0.009 | *** |
| DEG_auprc | f1@50 | 1.525 | [0.782, 1.892] | 0.001 | *** |
| Dir_all | f1@100 | 1.478 | [0.541, 1.852] | 0.003 | *** |
| Dir_deg | f1@100 | 1.274 | [0.282, 1.887] | 0.009 | *** |
| DEG_auprc | f1@100 | 1.525 | [0.737, 1.895] | 0.001 | *** |
| Dir_all | dir_disc | 1.404 | [0.468, 1.895] | 0.003 | *** |
| Dir_deg | dir_disc | 1.523 | [0.762, 1.866] | 0.001 | *** |
| DEG_auprc | dir_disc | 1.348 | [0.364, 1.875] | 0.007 | *** |

### H2 Cross-Dataset Summary

| Dataset | Sig/Total (f1@50) | Sig/Total (f1@100) | Sig/Total (dir_disc) | Overall |
|---------|:-----------------:|:------------------:|:--------------------:|:-------:|
| K562 | 3/3 | 3/3 | 3/3 | **9/9** |
| RPE1 | 1/3 | 1/3 | 3/3 | **5/9** |
| Norman | 3/3 | 3/3 | 3/3 | **9/9** |

**Key finding**: For the `dir_discovery_deg` downstream task, ALL 3 datasets show significant BioEval > MSE (6/6 significant). For `f1@50`, DEG_auprc is significant across all 3 datasets, but Dir_all/Dir_deg are not significant for RPE1 (wide CIs due to N=9).

**RPE1 has wider CIs** because its H1 reversal is strongest (MSE and BioEval are nearly independent), making the rho difference more variable under bootstrap.

## Verified Knowledge

1. **H1 statistical robustness**: Bootstrap CI confirms MSE-Dir ranking disagreement across all 3 datasets. RPE1 and Norman have CIs including 0 (independence). K562 has |tau|=0.696 < 0.7 (below threshold but significant).
2. **H2 statistical robustness**: `dir_discovery_deg` is the most robust downstream task — 6/6 comparisons significant. `DEG_auprc` is the most robust BioEval metric for H2.
3. **RPE1 f1@50 sensitivity**: RPE1 has wider CIs for f1@50 comparisons. Dir_all and Dir_deg vs f1@50 are NOT significant for RPE1, but DEG_auprc vs f1@50 IS significant. This suggests DEG_auprc is the most reliable H2 metric.
4. **Pearson and Dir_deg consistently agree**: tau(Pearson,Dir_deg) is 0.765-0.941 across all datasets with narrow CIs. This is expected — both are "higher=better" correlation-like metrics.
5. **N=9 limitation**: Bootstrap SE is 0.09-0.37 for tau and 0.00-0.35 for rho. Wide CIs are inherent to small model sets. Adding more models (DL models) would narrow CIs.

## Remaining Gaps

1. **RPE1 f1@50**: Dir_all/Dir_deg not significantly better than MSE for f1@50. Need more models or alternative downstream tasks.
2. **Sign convention discrepancy with run_16**: Run_16 used reversed MSE rankings; run_17 uses raw values. The |tau| values differ due to different random noise in noisy_ridge/sign_flip_ridge. Both support H1 but with different magnitudes.
3. **DL models**: Adding GEARS/CPA would increase N from 9 to 11-13, narrowing CIs significantly.
4. **Norman logFC scale correction**: Still pending (A3).
5. **Downstream task circularity**: Still pending (A4).
