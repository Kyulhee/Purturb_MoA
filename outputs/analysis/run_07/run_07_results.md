# Run 07 Results: Multi-Cell-Type FCR-ICM (RQ1 + RQ3 on Real Data)

**Date:** 2026-04-27 | **Status:** Completed

---

## Summary

FCR-ICM validated on real multi-cell-type data (Replogle 2022, K562+RPE1, 843 shared perturbations, 165K cells).

| Config | RQ1 (z_tx cross-CT corr) | RQ1 (cosine) | RQ3 (transfer R2) | RQ3 (transfer corr) |
|--------|--------------------------|--------------|--------------------|----------------------|
| FCR baseline (no ICM) | -0.349 | -0.319 | -0.297 | 0.411 |
| FCR + ICM | **0.348** | **0.376** | **0.924** | **0.967** |

## Key Findings

### 1. ICM is essential for cross-cell-type invariance on real data

Without ICM, z_tx is **negatively correlated** across cell types (corr=-0.349). The baseline encoder learns cell-type-specific perturbation effects — z_tx for the same perturbation points in opposite directions across K562 and RPE1.

With ICM, z_tx alignment flips to positive (corr=0.348). While lower than synthetic (0.971), this confirms ICM's real-data effectiveness.

### 2. RQ3: Zero-shot transfer succeeds dramatically with ICM

Transfer R2: -0.297 → 0.924 (+1.22). Transfer correlation: 0.411 → 0.967 (+0.56).

This is the most important result: using K562's z_tx to predict RPE1 expression achieves R2=0.92, nearly oracle-level performance.

### 3. RQ1 absolute values lower than synthetic

Synthetic: ICM corr=0.971. Real: ICM corr=0.348. Possible reasons:
- Real data has much higher dimensionality and noise
- 150 epochs may be insufficient for convergence with ICM
- icm_weight=10.0 may need tuning for real data
- Some perturbations may have genuinely different effects across cell types

### 4. Top perturbations achieve near-perfect alignment

With ICM, top perturbations reach corr 0.8-0.91:
- NOB1+ctrl: -0.715 → 0.910
- PTBP1+ctrl: -0.713 → 0.894
- TSR1+ctrl: -0.655 → 0.907

The per-perturbation variation suggests some genes are more cell-type-invariant than others.

### 5. RPE1 concatenate bug confirmed and fixed

The bug from run_04 (AnnData.concatenate(batch_key='cell_type') overwrites cell_type) was fixed by using batch_key='batch'. This was the root cause of the "RPE1 loading failure" — loading actually succeeded, but concatenation broke downstream logic.

## Experimental Details

- **Dataset:** Replogle 2022 K562 (162K cells) + RPE1 (161K cells)
- **Shared perturbations:** 843 (excl ctrl)
- **Shared genes:** 2,832 → 500 HVGs after preprocessing
- **Subsampling:** max 200 cells per perturbation per cell type
- **Training:** 150 epochs, Adam lr=1e-3, beta=0.5, icm_weight=10.0, z_dim=8
- **Evaluation:** 843 shared perturbations with ≥10 cells per cell type

## Implications

1. **H1 confirmed on real data**: ICM enables zero-shot cross-cell-type transfer (R2=0.92)
2. **RQ1/RQ3 real-data gap from stages/ resolved**: Both questions now have real-data evidence
3. **For the paper**: The transfer result (0.92 R2) is the headline finding — ICM unlocks zero-shot perturbation prediction across cell types
4. **Next steps**: Tune ICM weight/epochs to improve RQ1 absolute values; test on 4-cell-line dataset (Nadig/Replogle 2025)
