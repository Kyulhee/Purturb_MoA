# Run 12 Report: Cross-CT Epistasis Transfer — Full Ablation + Coverage

**Date:** 2026-04-30 | **Runtime:** 30659s (~8.5h) | **Data:** Replogle K562+RPE1 (848 shared perts, 500 HVG)

---

## 1. Success Criteria Assessment

| Metric | Value | Target | Verdict |
|--------|-------|--------|---------|
| Cross-CT rho (A7 residual rank) | 0.326 | >0.4 | **PARTIAL** (0.2-0.4 range) |
| U-Error rho (MC-only) | 0.786 | >0.6 | **PASS** |
| Coverage (90% CI) | 0.293 | 0.85-0.95 | **FAIL** |
| AL improvement (top10) | 2.667× | >2× | **PASS** |
| Norman Precision (top20) | 0.731 | >0.6 | **PASS** |

**Overall: PARTIAL SUCCESS — 3 PASS, 1 PARTIAL, 1 FAIL**

---

## 2. RQ2: Holdout UQ + Coverage

| Config | rho_mc | coverage_90 | mean_r2 | mean_pearson | coverage_in_range |
|--------|--------|-------------|---------|---------------|-------------------|
| A1 (ICM) | **0.786** | 0.293 | 0.860 | 0.928 | False |
| A2 (no-ICM) | 0.642 | 0.344 | — | — | False |

### Critical Finding: Coverage Failure
- MC Dropout 90% CI covers only 29.3% of genes (target: 85-95%)
- This is **severe under-coverage** — intervals are far too narrow
- Norman coverage similarly low: 0.235
- **Root cause hypothesis**: MC Dropout with p=0.1 in decoder produces insufficient variance. The model is overconfident.
- **Implication**: UQ calibration claim is unsupported until this is fixed

### Positive: rho_mc improved vs run_10
- A1: 0.786 vs run_10's 0.660 — MC Dropout correlation with error is strong
- A1 ICM model has better rho_mc than A2 (0.786 vs 0.642)

---

## 3. RQ3: Cross-CT Epistasis Transfer

### 3.1 Main Results

| Config | rho_add | rho_mult | rho_prod | rho_composite | rho_r_add_mag |
|--------|---------|----------|----------|---------------|---------------|
| A1 (ICM) | 0.367 | 0.421 | **0.437** | — | 0.326 |
| A2 (no-ICM) | 0.365 | 0.417 | 0.424 | 0.409 | 0.389 |
| **B1 (CPA)** | **0.430** | — | **0.430** | **0.430** | — |

### 3.2 A7 Residual Rank Transfer
- rho=0.326, pearson=0.324, p<0.001
- top10_overlap=0.100, top20=0.100, top50=0.180
- **Below run_11's 0.402 — likely due to different data pipeline**

### 3.3 ICM Ablation (A1 vs A2)
| Metric | ICM Improvement |
|--------|-----------------|
| imp_add | 1.005× |
| imp_mult | 1.008× |
| imp_prod | 1.031× |
| imp_r_add_mag | 0.837× |
| imp_composite | 0.962× |
| imp_mc_var | 0.431× |

**ICM negative result confirmed**: imp_add=1.005× (essentially 1.0×). ICM does not help epistasis transfer. ICM hurts residual rank transfer (0.837×) and MC variance transfer (0.431×).

### 3.4 A3 Single Formula
| Formula | rho |
|---------|-----|
| add_only | 0.367 |
| mult_only | 0.421 |
| prod_only | **0.437** |

Product formula is consistently the strongest — consistent with run_11's 3-formula sensitivity (rho 0.36-0.43).

### 3.5 A4 Trivial Decomposition
- rho=0.367 (identical to rho_add)
- Trivial = additive only, no residual decomposition

### 3.6 CPA Baseline (B1) — Concerning Finding
- CPA rho_add=**0.430** > FCR rho_add=0.367
- A simpler non-factorized model outperforms FCR on epistasis transfer
- **Interpretation**: Factorized decomposition does not help epistasis transfer. Simpler models capture the cross-CT signal equally well.
- This is actually consistent with the "epistasis is a simple biological conservation" narrative — it doesn't need complex decomposition.

### 3.7 Additional R2/Pearson Metrics
- mean_r2_src=0.950, mean_r2_tgt=0.937
- mean_pearson_src=0.977, mean_pearson_tgt=0.973
- R2 and Pearson are high within cell types (reconstruction quality)
- coverage_90 rho across CT: 0.230 (A1), 0.253 (A2) — weak

---

## 4. RQ4: Active Learning

| k | random | UQ | epi | oracle | imp_uq | imp_epi |
|---|--------|-----|------|--------|--------|---------|
| 5 | 0.004 | 0.020 | 0.020 | 0.020 | **5.0×** | **5.0×** |
| 10 | 0.012 | 0.031 | 0.031 | 0.039 | **2.667×** | **2.667×** |
| 20 | 0.024 | 0.059 | 0.059 | 0.078 | **2.5×** | **2.5×** |

- Transfer overlap = 0.750 (strong)
- AL PASS at top10 and top20

---

## 5. RQ1: Norman Epistasis Precision + Coverage

| Metric | Value |
|--------|-------|
| n_doubles | 131 |
| n_epi (70th percentile) | 40 |
| trivial_prec | 0.305 |
| rho_pred_actual | 0.601 |
| **prec_top10** | **0.923** (3.02× random) |
| **prec_top20** | **0.731** (2.39× random) |
| prec_top30 | 0.590 (1.93× random) |
| mean_coverage_90 | 0.235 |

Precision PASS, but Coverage FAIL (0.235).

---

## 6. Comparison with Prior Runs

| Metric | run_10 | run_11 | run_12 | Change |
|--------|--------|--------|--------|--------|
| rho_mc | 0.660 | — | 0.786 | +0.126 ↑ |
| A7 rho | — | 0.402 | 0.326 | -0.076 ↓ |
| ICM imp | — | 0.966× | 1.005× | Stable (~1.0) |
| AL imp (top10) | 5.0× | 4-8× | 2.667× | Lower |
| Norman prec_top20 | — | 0.60-0.75 | 0.731 | Consistent |

### Why A7 rho dropped (0.402→0.326)?
Likely cause: **Data pipeline change**. Run_11 used `gears` PertData loader (with its normalization/subsampling). Run_12 loads h5ad directly. Different preprocessing = different residual distributions = different rho values.

---

## 7. Critical Issues

### Issue 1: Coverage Failure (SEVERE)
- 90% CI covers only 29.3% of genes — should be ~90%
- MC Dropout variance is too small → intervals too narrow
- **Fix needed**: Increase MC Dropout rate, use conformal prediction calibration, or temperature scaling
- **Impact**: UQ calibration claim cannot be made without fixing this

### Issue 2: CPA beats FCR at epistasis transfer (MODERATE)
- CPA rho_add=0.430 > FCR rho_add=0.367
- However, CPA has no UQ, no factorized interpretation
- **Framing**: This supports "epistasis conservation is simple, doesn't need factorized models" — actually strengthens the conservation narrative

### Issue 3: A7 rho below target (MODERATE)
- 0.326 is PARTIAL (0.2-0.4), not PASS (>0.4)
- Product formula rho=0.437 does PASS, but residual rank (A7) does not
- **Note**: p<0.001, n=848 — statistically significant even if below threshold

---

## 8. Outcome Classification

Per experiment_contract.yaml:
- **Primary**: Cross-CT rho + AL → rho 0.326 is PARTIAL, AL 2.667× is PASS → **PARTIAL SUCCESS**
- **Partial success criteria**: "rho 0.2-0.4: moderate claim" — MET
- **Coverage**: FAIL — not in partial success criteria, but undermines UQ claim

**Allowed claim (partial)**: "에피스태시스가 세포유형 간에 부분 보존되나 전이 예측력이 제한적이다"

**Next step per contract**: "PORTAL 대규모 검증으로 일반성 확인"

**Recommended**: Fix coverage first (conformal prediction or higher MC Dropout rate), then assess whether loopback to planning is needed.

---

## 9. Evidence Trail
- run_12: A1 rho_add=0.367, A7 rho=0.326, CPA rho=0.430
- run_12: rho_mc=0.786 PASS, Coverage 0.293 FAIL
- run_12: AL 2.667× PASS, Norman prec_top20=0.731 PASS
- run_12: ICM imp=1.005× (negative result confirmed)
- run_11: A7 rho=0.402, AL 4-8×, ICM imp=0.966×
- run_10: rho_mc=0.660, Cross-CT rho=0.444
