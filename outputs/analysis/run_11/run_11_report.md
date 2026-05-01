# Run 11 Report: Cross-CT Epistasis Transfer + Full Ablation

**Date:** 2026-04-29 | **Elapsed:** 2887s (~48min)

---

## Executive Summary

Run 11 executed the full ablation design from Planning run_04. The **primary hypothesis fails**: ICM regularization does NOT improve cross-cell-type epistasis transfer (improvement ratio 0.966×, target >1.5×). However, epistasis IS cross-cell-type conserved (rho 0.36-0.43, all p<1e-26), and secondary metrics (AL, Norman Precision) pass their targets.

**Verdict: PRIMARY FAILURE — loopback to framing required per experiment_contract.yaml**

---

## 1. RQ2: Uncertainty Quantification

| Metric | Result | Target | Status |
|--------|--------|--------|--------|
| rho_combined | 0.524 | >0.6 | FAIL |
| rho_mc | 0.524 | — | MC-only dominates |
| rho_ztx | 0.147 | — | ICM signal weak for UQ |
| best weights | mc=1.0, ztx=0.0 | — | Combined collapses to MC-only |

**Assessment:** UQ rho=0.524 is below the 0.6 target. The ICM violation signal contributes nothing to UQ (weight=0.0 in optimal combination). This is consistent with the RQ3 finding that ICM alignment does not help for cross-CT tasks.

## 2. RQ3: Cross-CT Epistasis Transfer (PRIMARY)

### A1 (ICM model) vs A2 (no-ICM model)

| Signal | A1 (ICM) rho | A2 (no-ICM) rho | ICM Improvement |
|--------|-------------|-----------------|-----------------|
| add | 0.356 | 0.369 | **0.966×** |
| mult | 0.403 | 0.419 | **0.962×** |
| prod | 0.428 | 0.421 | 1.018× |
| composite | 0.363 | 0.370 | **0.983×** |
| icm_viol | 0.158 | 0.239 | 0.660× |
| r_add_mag | 0.402 | 0.413 | 0.974× |
| mc_var | 0.055 | 0.061 | 0.907× |

**Critical finding:** ICM regularization **slightly hurts** epistasis transfer on most metrics. The improvement ratio is <1.0 for all signals except prod (1.018×). ICM violation signal itself has poor cross-CT transfer (A1: 0.158 vs A2: 0.239), suggesting ICM alignment actually degrades the transferability of this signal.

### Top-k Overlap (Epistasis)

| k | A1 (ICM) | A2 (no-ICM) | Random | A1 Imp | A2 Imp |
|---|----------|-------------|--------|--------|--------|
| 10 | 0.100 | 0.100 | 0.012 | 8.4× | 8.4× |
| 20 | 0.100 | 0.200 | 0.024 | 4.2× | 8.4× |
| 50 | 0.320 | 0.320 | 0.059 | 5.4× | 5.4× |

A2 (no-ICM) actually has better top-20 overlap than A1 (ICM).

### A7/A8 Baselines vs A1

| Method | rho | Top-10 | Top-20 | Top-50 |
|--------|-----|--------|--------|--------|
| A1 (ICM transfer) | 0.356 | 0.100 | 0.100 | 0.320 |
| A7 (residual rank) | **0.402** | 0.100 | 0.050 | 0.220 |
| A8 (UQ rank) | 0.360 | 0.100 | 0.100 | 0.280 |

**A7 (simple residual rank transfer) outperforms A1 (ICM transfer) on rho (0.402 vs 0.356).** This is a devastating result for the ICM hypothesis.

### 3-Formula Cross-CT Sensitivity

| Formula | ICM rho | no-ICM rho | ICM > no-ICM? |
|---------|---------|-----------|---------------|
| add | 0.356 | 0.369 | No |
| mult | 0.403 | 0.419 | No |
| prod | 0.428 | 0.421 | Marginally |

Formula sensitivity is moderate (rho range 0.36-0.43). Product formula is most transferable.

## 3. RQ4: Active Learning

| k | Random | UQ | Epi | Oracle | imp_epi |
|---|--------|-----|-----|--------|---------|
| 5 | 0.000 | 0.016 | 0.016 | 0.020 | — (rand=0) |
| 10 | 0.008 | 0.032 | 0.032 | 0.040 | **4.0×** |
| 20 | 0.008 | 0.063 | 0.063 | 0.079 | **8.0×** |

- Transfer overlap: 0.800 (target >0.5) — **PASS**
- AL improvement: 4-8× (target >2×) — **PASS**

## 4. Norman Epistasis Precision

| Metric | Value | Target | Status |
|--------|-------|--------|--------|
| trivial_prec | 0.305 | — | Baseline |
| rho_pred_actual | 0.461 | — | Significant |
| prec_top10 | 0.750 | >0.6 | **PASS** |
| prec_top20 | 0.600 | >0.6 | **PASS** |
| prec_top30 | 0.500 | >0.6 | FAIL |

## 5. Overall Verdict vs Experiment Contract

### Primary Success Criteria

| Criterion | Result | Threshold | Verdict |
|-----------|--------|-----------|---------|
| Cross-CT rho | 0.356-0.428 | >0.4 | PARTIAL (add fails, prod passes) |
| ICM transfer improvement | 0.966× | >1.5× | **CRITICAL FAILURE** |

### Failure Criteria (from contract)

- "ICM 전이 < baseline: ICM이 에피스태시스 전이를 개선하지 못함" → **TRIGGERED**

### Contract-mandated next step: **Loopback to framing**

---

## 6. Interpretation

### What worked
1. **Epistasis IS cross-cell-type conserved** (rho 0.36-0.43, all p<1e-26). This is a genuine, non-trivial finding.
2. **AL works** (4-8× improvement, transfer overlap 0.8)
3. **Norman precision works** (0.60-0.75 vs trivial 0.30)
4. **Simple residual rank transfer (A7) is competitive** (rho=0.402)

### What failed and why
1. **ICM does NOT improve epistasis transfer.** The causal invariance constraint (MMD on z_tx across cell types) does not create a representation that transfers epistasis patterns better. Possible reasons:
   - ICM alignment forces z_tx to be cell-type-invariant, which may remove cell-type-specific epistasis signal
   - The MMD constraint optimizes for distribution matching, not for preserving perturbation-specific transferability
   - Epistasis may be inherently cell-type-specific enough that removing cell-type signal degrades transfer
2. **ICM violation signal has poor transferability** (rho=0.158 for ICM model vs 0.239 for no-ICM). Paradoxically, the ICM model's violation score is LESS transferable, suggesting the MMD constraint distorts the violation signal.
3. **A7 (simple residual magnitude rank) outperforms A1 (ICM transfer)**. The simplest baseline wins on rho.

### Implication for the primary claim
The claim "ICM 정렬이 교세포 에피스태시스 전이를 비인과 baseline 대비 유의미하게 개선한다" is **NOT supported by data**. The data supports a weaker claim: "에피스태시스가 세포유형 간에 부분 보존되나 ICM의 기여는 제한적" — which maps to the contract's "partial success" scenario, but even this is generous since ICM actively hurts on most metrics.

---

## 7. Comparison with Previous Runs

| Metric | run_10 | run_11 | Delta |
|--------|--------|--------|-------|
| Cross-CT rho | 0.444 | 0.356-0.428 | Similar range |
| RQ2 rho_combined | 0.660 | 0.524 | -0.136 (worse) |
| AL imp | 5.0× | 4.0-8.0× | Similar |
| Transfer overlap | 0.65 | 0.80 | +0.15 (better) |

The RQ2 degradation (0.660→0.524) may be due to different random seeds or the different gene count (500 vs previous).
