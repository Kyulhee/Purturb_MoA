# Run 10: Replogle Multi-Cell-Type Epistasis+UQ

**Date:** 2026-04-29 | **Dataset:** Replogle K562+RPE1 (843 shared perturbations, 500 HVGs)

---

## Key Result: RQ2 Target Achieved

**RQ2 rho_mc = 0.660 > 0.6 target** — run_09의 핵심 약점 해소

---

## Full Results

### RQ1: z_tx Cell-Type Invariance

| Config | Corr | Cosine |
|--------|------|--------|
| Baseline (no ICM) | -0.659 | -0.625 |
| ICM | 0.381 | 0.376 |

ICM이 z_tx를 교세포 양의 상관으로 반전. run_07과 일치.

### RQ3: Zero-Shot Cross-Cell-Type Transfer (K562→RPE1)

| Config | R2 | Corr | ztx_cc |
|--------|-----|------|--------|
| Baseline | -0.290 | 0.399 | -0.659 |
| ICM | **0.932** | **0.970** | 0.381 |

ICM으로 R2 -0.29→0.93 극적 개선. run_07과 일치.

### RQ2: Uncertainty Quantification

| Evaluation | Metric | Value | vs Target |
|-----------|--------|-------|-----------|
| **Holdout (non-circular)** | rho_mc | **0.660** | **>0.6 PASS** |
| Holdout | rho_ztx | 0.236 | — |
| Holdout | rho_combined | 0.660 | PASS |
| Holdout | Mean R2 | 0.868 | — |
| Composed pairs | rho_ood | -0.340 | FAIL |
| Composed pairs | rho_mc | 0.625 | — |
| Composed pairs | Mean residual mag | 0.250 | — |

**핵심 발견:**
- MC Dropout 분산이 단일 섭동 holdout에서 0.660 — run_09의 0.401에서 대폭 개선
- OOD distance는 composed pairs에서 음의 상관 (-0.340) — Replogle에서는 OOD가 유효하지 않음
- Best weight: MC=1.0, OOD/ztx=0.0 — MC Dropout만이 유효 UQ 신호

### Phase 5: Cross-CT Epistasis Consistency

| Metric | Value |
|--------|-------|
| Shared pairs | 200 |
| Cross-CT rho | 0.444 |
| Direction agreement | 0.590 |

K562와 RPE1 간 에피스태시스 랭킹이 유의미한 양의 상관(p<1e-10). 방향 일치 59% — 우연(50%)보다 유의미하나 강하지 않음.

### Phase 6: Active Learning

| Top-K | Random | OOD | Epi | imp_ood | imp_epi |
|-------|--------|-----|-----|---------|---------|
| 5 | 0.016 | 0.049 | 0.082 | 3.0x | **5.0x** |
| 10 | 0.066 | 0.066 | 0.164 | 1.0x | **2.5x** |
| 20 | 0.098 | 0.082 | 0.328 | 0.8x | **3.3x** |

- 에피스태시스 기반 AL: Top-5에서 5.0x 개선율 (run_09와 동일)
- OOD 기반 AL은 Top-10+에서 무효 (OOD 자체가 음의 상관이므로)
- **Transfer overlap = 0.65**: K562 Top-20 에피스태시스 중 65%가 RPE1에서도 강한 에피스태시스

---

## Run 09 vs Run 10 비교

| Metric | Run 09 (Norman) | Run 10 (Replogle) | 개선 |
|--------|----------------|-------------------|------|
| RQ2 rho (UQ) | 0.401 | **0.660** | +64% |
| RQ2 OOD rho | 0.385 | -0.340 | 악화 |
| RQ3 R2 (ICM) | 0.92 | 0.932 | 동일 |
| AL Top-5 epi | 5.0x | 5.0x | 동일 |
| Cross-CT | N/A | rho=0.444 | 신규 |
| AL Transfer | N/A | 0.65 overlap | 신규 |

---

## 검증된 지식

1. **RQ2 달성**: MC Dropout이 holdout 단일 섭동에서 rho=0.660. Norman(0.401) 대비 대폭 개선. Replogle이 더 큰 데이터(843 perts)에서 UQ 신호가 더 강함
2. **OOD 무효 (Replogle)**: OOD distance가 음의 상관. 이유 추정: 단일 섭동만 있어 composed z_tx가 모두 "OOD"인데, OOD≠에피스태시스. Norman은 double-KO가 일부 training에 있었으나 Replogle은 전무
3. **Cross-CT 에피스태시스 일관성**: K562↔RPE1 rho=0.444. 생물학적으로 합리적 — 공유 signaling pathway의 에피스태시스는 세포유형 간 부분 보존
4. **AL Transfer**: K562에서 식별한 Top-20 에피스태시스의 65%가 RPE1에서도 유효. ICM 정규화가 교세포 전이뿐 아니라 에피스태시스 전이에도 기여
5. **MC Dropout > OOD**: Replogle에서 MC Dropout이 유일한 유효 UQ 신호. OOD는 composed pairs에서 무효

## 남은 과제

- RQ2 OOD: Replogle에서 음의 상관. 설명 필요 (단일 섭동만 → composed z_tx가 모두 OOD)
- Cross-CT direction agreement 59%: 우연보다 높으나 강하지 않음. 세포유형 특이적 에피스태시스가 상당함
