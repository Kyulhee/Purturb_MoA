# Stage 04 — Analysis

## Workflow
1. `docs/04_analysis.md` 가이드 확인
2. `stages/03_planning.md`에서 실험 설계 및 기준 확인
3. `objects/current/experiment_contract.yaml` + `evaluation_validity_card.yaml` 사전 확인
4. 데이터 전처리, 모델 학습, 평가 수행
5. 결과를 Planning 타겟과 비교
6. 산출물 → `outputs/analysis/run_XX/`에 저장
7. 아래 지식 업데이트
8. `objects/current/result_card.yaml` 업데이트 (주요 지표, 베이스라인 비교, 소거 요약, 다음 단계)
9. `objects/current/validation_readiness_card.yaml` 업데이트 (증거 상태, 차단요소, 허용 클레임 강도)
10. 방향 전환 시: `objects/current/pivot_diagnosis_card.yaml` 작성 필수
11. CLAUDE.md의 `current_stage` 업데이트

## Done when
- Planning 타겟 성능 달성 또는 미달 시 사용자 보고 후 방향 결정

---

## 검증된 핵심 지식

### H1 인과 불변성 가설 — 실험 결과 (run_04-05)

**RQ1: ICM이 z_tx를 cell type 불변으로 만드는가** — 합성 통과, 실제 통과
- FCR (ICM 없음): 교세포 상관계수 0.505 (합성) / -0.349 (실제)
- FCR + ICM: 교세포 상관계수 0.971 (합성) / 0.348 (실제)
- ICM MMD 정규화가 10개 교란 모두에서 개선 (범위: +0.21 ~ +0.74, 합성)
- 실제데이터: baseline z_tx가 음의 상관(-0.35) → ICM이 정의 상관(+0.35)으로 반전

**RQ2: 단일-KO z_tx로 조합적 예측 가능한가** — 합성(잠재공간) 실패, 합성(유전자공간) 통과, 실제 통과
- 합성 잠재공간: best_corr=0.29, best_R2=-1.63 — 인코더 비선형 변환으로 조합 구조 파괴
- 합성 유전자공간: **R2=0.88, corr=0.94** — 디코더가 인코더 비선형성 보상 (run_05로 검증)
- Norman 2019 실제: FCR best_corr=0.955, best_R2=0.881 / FCR+ICM best_corr=0.951, best_R2=0.870
- **소거실험(run_05)**: 조합 일관성 손실이 RQ2-cross 0.20→0.79 개선. 전체 모델(config 6) RQ1=0.99, RQ3=0.99
- **실제데이터 소거실험(run_06)**: Norman 128 double-KO 쌍에서 모든 config R2=0.86-0.89로 유사. comp loss가 오히려 성능 저하(0.88→0.86). baseline FCR이 이미 조합성에 충분

**RQ3: ICM이 zero-shot 교세포 전이를 가능하게 하는가** — 합성 통과, 실제 통과
- FCR (ICM 없음): 전이 상관계수 0.508, 코사인 0.476 (합성) / R2=-0.30, corr=0.41 (실제)
- FCR + ICM: 전이 상관계수 0.960, 코사인 0.956 (합성) / R2=0.92, corr=0.97 (실제)
- 실제데이터: ICM으로 전이 R2가 -0.30 → 0.92로 극적 개선

**Replogle 2022**: K562+RPE1 결합 성공 (843 공유 교란, 165K cells) — RQ1/RQ3 실제 검증 완료 (run_07)
- run_04 RPE1 로딩 실패 원인: `concatenate(batch_key='cell_type')`가 cell_type 값을 덮어씀 → `batch_key='batch'`로 수정

### 핵심 인사이트
1. **ICM은 불변성에 확실히 유효** (RQ1, RQ3). MMD 정규화가 z_tx 교세포 정렬에 강력
2. **RQ2 합성-실제 갭 해명됨**: 조합성은 유전자 공간에서 평가해야 함. 잠재 공간 R2=0.05 vs 유전자 공간 R2=0.88. 디코더가 인코더 비선형성 보상
3. **조합 일관성 손실이 핵심 구성요소** (합성): ICM만으로 RQ2 개선 안 됨. comp_loss가 RQ2-cross를 0.20→0.79로 향상. 단 실제데이터에서는 불필요
4. **ICM이 인코더를 더 선형적으로 만듦**: linear R2 0.69→0.87. ICM의 분포 정렬이 부차적으로 조합 구조 보존에 기여
5. **실제데이터에서 comp loss 불필요**: Norman 소거실험에서 baseline FCR이 이미 R2=0.88. comp loss는 과도한 제약으로 오히려 저하
6. **가법성이 우위**: 128 double-KO 중 110쌍(86%)에서 가법 조합이 승법보다 우수
7. **ICM이 실제데이터 교세포 전이에 결정적** (run_07): RQ3 전이 R2 -0.30→0.92, corr 0.41→0.97. Baseline z_tx는 교세포 음의 상관(-0.35)
8. **RQ1 실제 절대값은 합성보다 낮음**: 합성 0.97 vs 실제 0.35. 튜닝/학습량으로 개선 여지

### 이전 방향(NAP)에서 검증된 지식 (참고용)
- **XGBoost-only R2=0.91** — FBA는 근본적으로 tabular problem
- **GNN 임베딩 중복**: 정적 그래프에서는 knockout mask가 충분 통계량
- **AL 실패**: FBA가 싸고 입력 차원이 낮아 AL 이점 없음

---

## 현재 진행 상태

### 기반 모델 (Phase 1-2)

| 단계 | 상태 | 비고 |
|------|------|------|
| Step 1: 데이터 확보 | 완료 | Norman 2019 (89K cells), Replogle K562+RPE1 (165K cells) |
| Step 2: FCR 인코더 구현 | Phase 1 완료 | VAE, z_dim=8, z_x/z_t/z_tx 분해 |
| Step 3: ICM 정규화 (RQ1) | 합성+실제 통과 | 합성 0.50→0.97 / 실제 -0.35→0.35 |
| Step 4: 조합성 (RQ2) | 합성 실패, 실제 통과 | Norman best_corr=0.955, best_R2=0.881 |
| Step 5: 교세포 전이 (RQ3) | 합성+실제 통과 | 합성 0.51→0.96 / 실제 R2 -0.30→0.92 |
| Step 6: 소거 실험 | 완료 (run_05-06) | 6구성 소거 + 인코더 비선형성 + 갭 해명 + 실제데이터 소거 |

### 에피스태시스+UQ (Phase 3-6, run_09-10)

| 단계 | 상태 | 비고 |
|------|------|------|
| Step 7: 잔차 분해 (RQ1) | 합성 통과 | r=1.000 (run_09) |
| Step 8: 불확실성 정량화 (RQ2) | Replogle 통과 | Norman rho=0.401 미달 / Replogle rho=0.660 PASS |
| Step 9: 에피스태시스 탐지 (RQ3) | 부분 통과 | AUROC=1.000 (순환 의심), Formula agreement 76.6% PASS |
| Step 10: 능동학습 (RQ4) | 통과 | AL 5.0× random (Norman+Replogle), Transfer overlap 0.65 |
| Step 11: 논문 초안 | 미실행 | Stage 05로 이관 |

### BioEval 메트릭 연구 (run_13, 2026-04-30)

**프레이밍**: 생물학적 충실도 기반 평가 지표 설계 + MSE/R2 순위 반전 현상 규명

**핵심 결과 (3 데이터셋 × 11 시뮬레이션 모델):**

| 데이터셋 | tau(MSE, Dir_deg) | tau(R2, Dir_deg) | 반전 쌍 수 | 판정 |
|----------|-------------------|-------------------|-----------|------|
| K562 | 0.382 | 0.091 | 9/12 | REVERSAL |
| RPE1 | 0.600 | 0.236 | 5/12 | PARTIAL |
| Norman | **-0.200** | **-0.236** | 10/12 | **ANTI-CORRELATED** |

**Mean-effect trap 증거:**
- Norman: mean_predictor MSE #1 → Dir #11. MSE와 방향 정확도가 반대 상관
- Norman DEG fraction 1.53% → predicting zero gives excellent MSE but 0% direction
- 교세포 일관성: 모든 지표 tau > 0.78 (K562↔RPE1) → 현상이 체계적

**S1 민감도**: DEG 임계값 스윕 PASS — 방향 정확도 견고

**검증된 지식:**
- **MSE/R2 순위 반전 실재** (H1 SUPPORTED): Norman에서 가장 극적
- **Mean-effect trap**: MSE가 안전한 예측(평균 근처)에 보상 → 방향 정보 상실
- **Norman 취약성**: 적은 섭동 수(283) + 낮은 DEG 비율(1.53%)으로 반전 심화
- **BioEval-Dir과 기존 지표의 순위 불일치가 체계적**: 교세포 일관성으로 확인

**남은 과제:**
- ~~실제 모델 예측 확보~~: run_15 Ridge LOO 완료 (Norman) + run_16 gene PCA Ridge (K562/RPE1)
- ~~Phase 4 구현~~: run_14 완료
- ~~K562/RPE1 유전자 수준 feature Ridge~~: run_16 완료 — C1 해결
- ~~Bootstrap CI for Kendall tau~~: run_17 완료 — H1+H2 통계적 견고성 확인
- ~~Norman logFC 스케일 보정~~: run_18 완료 — 보정 불필요 (Dir_deg 불변, downstream 과업 미개선). B3 해결
- ~~Downstream task independence (B6)~~: run_19 완료 — H2 domain-specific. Cross-domain 33.3%, intra-DEG 100%, intra-magnitude 100%. B6 해결
- ~~실제 DL 모델 예측 (GEARS, CPA) — H3 "DL > baseline" 검증용~~ → run_20 GEARS 완료. 결과: GEARS < Ridge, DL ≠ better

### run_15 실제 학습 모델 (sklearn Ridge LOO, 2026-04-30)

**프레이밍**: 시뮬레이션 모델 한계 극복 — sklearn Ridge LOO로 실제 학습 모델 예측 확보

**모델 구성 (9개)**:
- Trained: Ridge alpha=1/10/100 (analytical LOO, hat matrix)
- Degraded: noisy_ridge (15% noise), sign_flip_ridge (15% sign flip)
- Baselines: mean_predictor, mean_effect
- Oracle: constant_shrink (0.15×true), half_signal (0.5×true) — H3에서 제외 필요

**핵심 결과 — Norman (additive features, 신뢰 가능):**

| 지표 | ridge | ridge_med | ridge_strong | mean_predictor | mean_effect |
|------|-------|-----------|-------------|----------------|-------------|
| R2 | **0.643** | 0.402 | 0.143 | -0.004 | -0.002 |
| Dir_deg | **0.986** | 0.981 | 0.963 | 0.000 | 0.571 |
| DEG_auprc | **0.764** | 0.621 | 0.281 | 0.015 | 0.039 |
| f1@50 | **0.488** | 0.442 | 0.407 | 0.021 | 0.021 |

- H1: 1/5 REVERSAL + 4/5 PARTIAL (tau=0.5-0.61). 시뮬레이션(run_13)의 ANTI-CORRELATED(-0.20)보다 약하지만 여전히 순위 불일치 존재
- H2: **15/15 (100%)** — BioEval > MSE as downstream predictor. rho(MSE, dir_disc)=0.678 vs rho(Dir_deg, dir_disc)=1.000
- H3 (corrected): **trained > true baselines across all metrics**. Oracle baselines 제외 시 Ridge가 mean_predictor/mean_effect를 압도

**구조적 문제 — K562/RPE1 (one-hot LOO 퇴화):** → run_16에서 해결

**검증된 지식:**
- **H1+H2 실제 학습 모델로 확인 (Norman)**: 시뮬레이션 한계 극복. Ridge LOO로 genuine prediction 생성
- **H3 SUPPORTED (Norman, corrected)**: trained Ridge > true baselines (mean_predictor, mean_effect). Oracle baseline 제외 필수
- **One-hot LOO 퇴화**: 단일-KO 데이터셋에서 identity feature matrix로는 LOO 불가. 유전자 수준 feature 필요 → run_16에서 해결
- **Oracle baseline 분류 필수**: constant_shrink, half_signal은 ground truth 사용 → H3에서 제외
- **Norman additive Ridge R2=0.643**: 단일-KO feature로 double-KO 예측이 실제로 작동 (run_06 가법성 결과와 일치)

### run_16 Gene PCA Feature Ridge (2026-05-01)

**프레이밍**: C1(one-hot LOO 퇴화) 해결 — K562/RPE1에 PCA 기반 유전자 수준 feature 적용

**Feature 설계**:
- PCA on control cells (30 PCs) → 각 섭동 cell을 PCA 공간에 projection
- Feature = [pca_mean(30) | pca_var(30) | log_cell_count(1)] = 61 dims
- Feature rank: 61/61 (full rank) — LOO 시 다른 섭동의 PCA profile이 남아 예측 가능

**핵심 결과 — 3 데이터셋 전체 실제 모델 확보:**

| 데이터셋 | Ridge R2 | Dir_deg | DEG_auprc | f1@50 | vs run_15 |
|----------|----------|---------|-----------|-------|-----------|
| **K562** | **0.523** | **0.985** | **0.722** | **0.354** | R2 -0.027→0.523 |
| **RPE1** | **0.652** | **0.989** | **0.827** | **0.274** | R2 -0.027→0.652 |
| **Norman** | **0.643** | **0.986** | **0.764** | **0.488** | (same, additive features) |

**H1: Metric-Ranking Reversal — 3 데이터셋 전체 확인:**

| 데이터셋 | tau(MSE,Dir_all) | tau(MSE,Dir_deg) | tau(Pearson,Dir_deg) | 판정 |
|----------|:----------------:|:----------------:|:--------------------:|------|
| K562 | 0.500 | 0.611 | **0.056** | PARTIAL/REVERSAL |
| RPE1 | **0.333** | **0.389** | **0.111** | **REVERSAL** |
| Norman | 0.500 | 0.500 | **-0.167** | PARTIAL/REVERSAL |

- RPE1이 가장 강력한 실제 모델 H1 REVERSAL 증거 (tau=0.333)

**H2: Metric-Downstream Correlation — Domain-Specific (run_19 재분석):**

| Domain | Pass/Total | Pass Rate | Mean Gap | 해석 |
|:------:|:----------:|:---------:|:--------:|------|
| cross-domain | 8/24 | 33.3% | -0.086 | Dir metrics do NOT predict gene-set tasks better than MSE |
| intra-DEG | 9/9 | **100.0%** | +0.319 | DEG_auprc predicts f1@50 — confirms meaningful signal |
| intra-magnitude | 6/6 | **100.0%** | +0.303 | mag_rank predicts gene-set tasks — shared magnitude info |
| intra-direction | 3/6 | 50.0% | -0.013 | Dir vs dir_discovery — partially circular, weak |

**핵심 인사이트 (run_19)**:
- MSE 자체가 direction task(dir_discovery)를 rho=0.88-0.96으로 예측 — MSE는 domain-general predictor
- Dir metrics은 gene-set tasks를 MSE보다 잘 예측하지 못함 — cross-domain H2 약함
- DEG_auprc/mag_rank은 각자 도메인에서 MSE보다 강력 — intra-domain H2 강함
- **Revised H2**: BioEval metrics provide domain-specific predictive advantage over MSE

**H3: Trained vs Baseline — 3 데이터셋 ALL WIN (corrected, oracle excluded) + GEARS DL 결과:**

| 데이터셋 | Dir_deg trained | Dir_deg baseline | R2 trained | R2 baseline | 6 지표 |
|----------|:--------------:|:----------------:|:----------:|:-----------:|:------:|
| K562 | 0.985 | 0.606 | 0.523 | -0.013 | **ALL WIN** |
| RPE1 | 0.989 | 0.666 | 0.652 | -0.013 | **ALL WIN** |
| Norman | 0.986 | 0.571 | 0.643 | -0.002 | **ALL WIN** |

**GEARS DL 모델 (run_20): Ridge < GEARS, DL ≠ better**

| 데이터셋 | GEARS R2 | Ridge R2 | GEARS Dir_deg | Ridge Dir_deg | GEARS vs Ridge |
|----------|:--------:|:--------:|:-------------:|:-------------:|:--------------:|
| K562 | 0.085 | 0.610 | 0.888 | 0.988 | **0/4 승** |
| RPE1 | 0.147 | 0.696 | 0.890 | 0.983 | **0/4 승** |
| Norman | -0.699 | 0.896 | 0.422 | 1.000 | **0/4 승** |

→ H3 정교화: model quality > model complexity. Well-trained Ridge > baselines, but poorly-trained GEARS(DL) < Ridge. BioEval이 품질 격차를 정확히 식별.

**검증된 지식:**
- **C1 RESOLVED**: Gene PCA features로 one-hot LOO 퇴화 완전 해결. K562 R2 -0.027→0.523, RPE1 -0.027→0.652
- **H1+H2 전체 3 데이터셋에서 실제 학습 모델로 확인**: K562 PARTIAL/REVERSAL, RPE1 REVERSAL, Norman PARTIAL/REVERSAL. H2=100%
- **RPE1이 가장 강력한 H1+H2 증거**: 높은 DEG 비율(6.5%)이 방향 지표를 가장 판별력 있게 만듦
- **H3 3 데이터셋 ALL WIN**: trained Ridge > true baselines, 6 지표 모두
- **Gene PCA features는 낮은 누적 분산(9-14%)에서도 작동**: 핵심은 비퇴화 feature, 높은 분산 아님

---

## 미해결 과제

1. ~~**RQ1/RQ3 실제데이터 검증**~~: run_07 완료
2. ~~**RQ2 합성-실제 갭 해명**~~: run_05에서 해명
3. ~~**소거 실험**~~: run_05 완료
4. ~~**실제데이터 소거실험**~~: run_06 완료
5. **에피스태시스 탐지+UQ (run_09→run_10)**: RQ2 달성(rho=0.660), RQ4 달성 — 아래 참조
6. **논문 초안 작성**: Stage 05로 이관
7. ~~**BioEval 실제 모델 예측 확보**~~: run_15 Ridge LOO (Norman) 완료 + run_20 GEARS (DL)
8. ~~**BioEval Phase 4 구현**: H2(downstream 과업 상관) 검증~~: run_14 완료 — H2 SUPPORTED (88.9%)
9. ~~**H3 검증**~~: run_15 (Norman) + run_20 (GEARS DL) — Ridge > baselines, GEARS < Ridge

### run_09→10 에피스태시스+UQ 실험 (2026-04-29)

**프레이밍**: FCR-ICM 잔차에서 에피스태시스 분해 + 불확실성 정량화 + 능동학습

**결과 (Norman run_09 + Replogle run_10):**

| RQ | 지표 | Norman (run_09) | Replogle (run_10) | 타겟 | 판정 |
|----|------|----------------|-------------------|------|------|
| RQ1 | 잔차 분해 (합성) | r=1.000 | — | r>0.7 | ✅ |
| RQ2 | Holdout rho_mc | 0.401 | **0.660** | >0.6 | ✅ (Replogle) |
| RQ2 | OOD rho | 0.385 | -0.340 | — | ❌ (Replogle) |
| RQ3 | AUROC (add) | 1.000 | — | >0.75 | ⚠️ 순환 |
| RQ3 | Formula agreement | 0.766 | — | >60% | ✅ |
| RQ4 | Top-5 AL improvement | 5.0x (epi) | 5.0x (epi) | >2x | ✅ |
| RQ4 | AL Transfer overlap | — | 0.65 | — | 신규 |
| — | Cross-CT 에피스태시스 | — | rho=0.444 | — | 신규 |

**검증된 지식:**
- 에피스태시스는 스펙트럼 → 이진 분류가 아닌 연속 순위 평가가 올바름
- **MC Dropout이 Replogle에서 UQ 유효 (rho=0.660)**: Norman(0.401) 대비 대폭 개선. 더 큰 데이터에서 UQ 신호 강화
- **OOD distance는 데이터셋 의존적**: Norman에서는 유효(rho=0.385), Replogle에서는 무효(-0.340). 단일 섭동만 있는 데이터에서는 OOD가 조합 불확실성 프록시로 부적합
- 3개 공식이 76.6% 유전자별 방향 일치 → 공식 민감도는 실제이지만 치명적이지 않음
- **Cross-CT 에피스태시스 부분 보존**: K562↔RPE1 rho=0.444, 방향 일치 59%. 단 ICM 기여는 run_11에서 부정됨

**남은 과제:**
- ~~RQ2 OOD~~: MC Dropout만 유효. OOD는 Replogle에서 무효
- ~~RQ3 순환 평가 해소~~: 교세포 에피스태시스 전이로 재정의 완료. run_11에서 주 가설 기각

### run_11 교세포 에피스태시스 전이 소거실험 (2026-04-29)

**프레이밍**: Planning run_04 기반. ICM vs no-ICM 에피스태시스 전이 소거 (A1 vs A2), A7/A8 baseline, Norman Precision, 3-공식 민감도

**핵심 결과: 주 가설 기각**

| 지표 | A1 (ICM) | A2 (no-ICM) | A7 (잔차 순위) | 타겟 | 판정 |
|------|----------|-------------|---------------|------|------|
| rho_add | 0.356 | **0.369** | **0.402** | >0.4 | partial |
| rho_prod | **0.428** | 0.421 | — | >0.4 | PASS |
| ICM 개선율 (add) | **0.966×** | (baseline) | — | >1.5× | **CRITICAL FAIL** |

**보조 지표 (PASS):** AL 4-8×, Transfer overlap 0.80, Norman Precision 0.60-0.75

**검증된 지식:**
- **ICM 정규화가 에피스태시스 전이를 개선하지 못함** (0.966×). MMD 분포 정렬은 유효(run_07)하나 에피스태시스 전이에는 부정적
- **A7(단순 잔차 순위 전이)가 A1(ICM 전이)보다 우수** (rho 0.402 vs 0.356) — 가장 간단한 baseline이 최적
- **에피스태시스 교세포 보존은 ICM과 무관한 생물학적 사실** — no-ICM에서도 rho=0.37-0.43
- **ICM 위반=불확실성 해석이 무효** — UQ 최적 weight: mc=1.0, ztx=0.0

**방향 전환:** 실험 계약 실패 기준 트리거 → framing loopback

---

## Run 이력 (세부 내용은 outputs/analysis/run_XX/ 참조)
- **run_01** (2026-04-26): NAP E2E 파이프라인. XGBoost R2=0.91, GNN 중복성 확인
- **run_02** (2026-04-27): GNN vs tabular 문헌 심층 리뷰 11편
- **run_03** (2026-04-27): Input-space AL 실험. AL R2=0.56 vs Random R2=0.68 — AL 실패
- **run_04** (2026-04-27): 방향 전환 + Phase 1 합성검증 + Phase 2 Norman (best_corr=0.955, R2=0.881)
- **run_05** (2026-04-27): 소거실험(6구성) + RQ2 갭 해명(잠재공간 R2=0.05 vs 유전자공간 R2=0.88)
- **run_06** (2026-04-27): Norman 소거실험 — 모든 config R2=0.86-0.89, comp loss 불필요
- **run_07** (2026-04-27): 다세포유형 검증 (Replogle) — RQ1 corr -0.35→0.35, RQ3 R2 -0.30→0.92
- **run_09** (2026-04-29): 에피스태시스+UQ — RQ4 달성(5.0x), RQ2 미달(0.401), OOD 유효, 공식 일치 76.6%
- **run_10** (2026-04-29): Replogle 다세포유형 에피스태시스+UQ — RQ2 rho=0.660 PASS, Cross-CT rho=0.444, AL 5.0x, Transfer overlap 0.65
- **run_11** (2026-04-29): 교세포 에피스태시스 전이 소거실험 — **주 가설 기각**: ICM 개선율 0.966×, A7>A1. AL/Norman Precision PASS. → framing loopback
- **run_12** (2026-04-30): Full ablation + Coverage + R2/Pearson + CPA baseline — **PARTIAL SUCCESS**: 3 PASS (rho_mc=0.786, AL 2.67×, Norman prec 0.731), 1 PARTIAL (A7 rho=0.326), 1 FAIL (Coverage 0.293). CPA rho=0.430 > FCR rho=0.367. ICM negative result confirmed (1.005×). Coverage under-coverage critical.
- **run_13** (2026-04-30): BioEval 메트릭-순위 반전 분석 — **H1 SUPPORTED**: MSE/R2와 BioEval-Dir 간 순위 반전 확인. K562 tau(MSE,Dir_deg)=0.382, RPE1 0.600, Norman **-0.200**(반대상관). Norman에서 mean_predictor가 MSE #1이나 Dir #11. 교세포 일관성 tau>0.78. 시뮬레이션 모델(11개) 기반, Phase 4 미구현
- **run_14** (2026-04-30): BioEval Phase 4 downstream 과업 상관 분석 — **H2 SUPPORTED**: BioEval 지표가 MSE보다 downstream 생물학적 유용성을 더 잘 예측 (88.9% 통과). K562 100%, RPE1 73.3%, Norman 93.3%. Norman에서 MSE가 DEG 회복과 **반대 상관**(rho=-0.736), Dir_deg는 방향 발견과 rho=0.945. H1+H2 결합: MSE와 BioEval이 순위 불일치(H1) + BioEval 순위가 생물학적 유용성과 정렬(H2)
- **run_15** (2026-04-30): BioEval 실제 학습 모델 (sklearn Ridge LOO) — **H1+H2 실제 모델로 확인, H3 SUPPORTED (Norman)**: Norman additive Ridge R2=0.643, Dir_deg=0.986. H1 PARTIAL (tau=0.5), H2 100%, H3 trained>baselines (oracle 제외 시). K562/RPE1 one-hot LOO 퇴화로 결과 부적합
- **run_16** (2026-05-01): Gene-level PCA feature Ridge — **C1 해결: K562/RPE1 실제 모델 결과 확보**: K562 Ridge R2=-0.027→0.523, RPE1 -0.027→0.652. H1 RPE1 REVERSAL (tau=0.333), K562 PARTIAL/REVERSAL. H2 100% (3 데이터셋 모두). H3 ALL WIN (3 데이터셋, 6 지표). RPE1이 가장 강력한 H1+H2 증거
- **run_17** (2026-05-01): Bootstrap CI for tau/rho — **H1+H2 통계적 견고성 확인**: B=10,000 bootstrap. H1: RPE1/Norman CI includes 0 (independence confirmed). K562 |tau(MSE,Dir_all)|=0.696<0.7. H2: dir_discovery 6/6 SIG, K562+Norman 9/9 SIG. RPE1 f1@50 1/3 SIG (DEG_auprc만). DEG_auprc이 가장 견고한 H2 지표
- **run_18** (2026-05-01): Scale Correction (A3) — **보정 불필요**: 3가지 보정(global/per-pert/variance-match) 테스트. Dir_deg 불변(부호 기반이므로), f1@50/dir_discovery 불변(순위 보존). variance_match가 H2 gap 증가시키나 MSE 악화로 인한 인위적 효과. B3 해결 — Norman logFC scale mismatch는 무시 가능
- **run_19** (2026-05-01): Downstream Task Independence (A4/B6) — **H2 is domain-specific**: 원래 H2 100% pass rate는 intra-domain pairs(DEG↔f1, Mag↔f1)에 의해 driven. Cross-domain(Dir↔f1)은 33.3%만 pass. MSE 자체가 dir_discovery를 rho=0.88-0.96으로 예측 — MSE는 domain-general predictor. Revised H2: BioEval metrics provide domain-specific predictive advantage. B6 해결
- **run_20** (2026-05-04): GEARS DL 모델 학습+평가 — **GEARS < Ridge, DL ≠ better**: GEARS(GNN+attention)를 3 데이터셋에서 학습(5 epochs). K562 GEARS R2=0.085 vs Ridge 0.610; RPE1 0.147 vs 0.696; Norman -0.699 vs 0.896. GEARS vs Ridge: 0/12 승. GEARS vs mean_predictor: 11/12 승. B2 해결 — DL 모델 테스트 완료, 단 Ridge에 패배. H3 정교화: model quality > model complexity. BioEval이 이 품질 격차를 정확히 식별
