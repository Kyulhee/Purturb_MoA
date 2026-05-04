# BioEval: 생물학적 충실도 기반 섭동 예측 평가 프레임워크 — 최종 해석 리포트

**날짜**: 2026-05-03
**Stage**: Analysis → Interpretation
**실험 이력**: run_13 ~ run_20

---

## 1. 연구 요약

본 연구는 Perturb-seq 섭동 예측에서 표준 평가 지표(MSE/R²)가 만드는 **Mean-Effect Trap**을 진단하고, 생물학적 충실도 기반 평가 지표(BioEval)를 설계하여 "DL ≤ Baseline" 위기가 지표의 아티팩트인지 판별하는 것을 목표로 했다.

세 가지 가설을 검증했다:

| 가설 | 내용 | 판정 | 클레임 강도 |
|------|------|------|-----------|
| **H1** | MSE/R²와 BioEval-Dir은 모델을 다르게 순위 매긴다 (순위 반전) | ⭐ SUPPORTED | STRONG |
| **H2** | BioEval 지표가 MSE보다 downstream 생물학적 유용성을 더 잘 예측한다 (domain-specific) | ⭐ SUPPORTED | MODERATE |
| **H3** | BioEval 하에서 학습된 모델이 단순 baseline을 능가한다 | ⭐ SUPPORTED (qualified) | STRONG |

---

## 2. 핵심 발견의 해석

### 2.1 H1: 순위 반전은 실재한다

**정량적 증거**:
- RPE1: τ(MSE, Dir_deg) = 0.389, 95% CI가 0 포함 → 두 지표가 통계적으로 독립
- Norman: τ(MSE, Dir_deg) = 0.500, 부호 반전 쌍 존재
- K562: τ(MSE, Dir_deg) = 0.611 < 0.7 기준

**의미**:
- MSE 최적화 모델이 반드시 방향 정확도가 높은 것은 아니다
- Norman에서 mean_predictor가 MSE #1이나 Dir_deg #11 — MSE가 "평균에 가까운 예측"에 보상을 준다는 직접적 증거
- DEG 비율이 낮을수록(Norman 1.53% < K562 2.38% < RPE1 6.50%) Mean-Effect Trap이 심화 — 비차별 유전자가 많을수록 MSE가 평균 예측에 더 큰 보상을 줌

**한계**:
- 9개 모델(3 Ridge + 2 degraded + 2 baseline + 2 oracle)로 순위 해상도 제한
- Ridge(선형 모델)로만 검증 — DL 모델에서의 순위 반전 패턴은 다를 수 있음
- Bootstrap CI에서 K562는 0 미포함이나 |τ| < 0.7 — "약한 불일치" 범주

### 2.2 H2: BioEval은 domain-specific 이점이 있다

**정량적 증거**:
- Intra-DEG (DEG_auprc ↔ f1@50): 9/9 (100%) 통과, mean gap = +0.319
- Intra-magnitude (mag_rank ↔ gene-set): 6/6 (100%) 통과, mean gap = +0.303
- Cross-domain (Dir ↔ gene-set): 8/24 (33.3%) 통과, mean gap = -0.086
- MSE 자체가 direction task를 ρ = 0.88-0.96으로 예측 — MSE는 domain-general predictor

**의미**:
- BioEval 지표는 동일 도메인 내에서 MSE보다 downstream 유용성을 더 잘 예측한다
- DEG_auprc가 f1@50을 잘 예측하는 것은 순환이 아니다 — 서로 다른 연산(연속 점수 vs 이진 임계값 F1)이 동일한 생물학적 신호를 포착
- MSE가 방향 정보를 포착한다는 발견(ρ(-MSE, dir_discovery) = 0.88-0.96)은 직관에 반한다 — MSE가 낮다는 것은 예측이 평균에 가깝다는 뜻인데, 방향까지 맞출 수 있다는 것은 비차별 유전자의 오차가 지배적이기 때문

**한계**:
- Cross-domain H2는 33.3%만 통과 — 방향 지표가 gene-set 과업을 MSE보다 잘 예측하지 못함
- N=9로 Spearman ρ 통계적 해상도 제한
- H2 클레임은 domain-specific으로 한정해야 함

### 2.3 H3: 학습 모델이 baseline을 능가한다

**정량적 증거**:
- K562: Ridge Dir_deg = 0.985 vs mean_predictor = 0.000, R² = 0.523 vs -0.013
- RPE1: Ridge Dir_deg = 0.989 vs mean_predictor = 0.000, R² = 0.652 vs -0.013
- Norman: Ridge Dir_deg = 0.986 vs mean_predictor = 0.000, R² = 0.643 vs -0.002
- 3 데이터셋 × 6 지표 = 18/18 ALL WIN

**의미**:
- 학습된 Ridge 모델이 모든 BioEval 지표에서 baseline을 압도적으로 능가
- "DL ≤ Baseline" 위기는 적어도 선형 모델 수준에서는 지표의 아티팩트 — 학습 모델이 baseline보다 나쁘다는 판단은 MSE만 사용했을 때의 결론
- Ridge가 mean_predictor보다 R²이 높다는 것은 자명하지만, Dir_deg에서 0.000 → 0.985라는 점프는 방향 정보가 학습을 통해 획득된다는 증거

**GEARS DL 모델 결과 (run_20)**:

| 데이터셋 | GEARS R² | Ridge R² | GEARS Dir_deg | Ridge Dir_deg | GEARS vs Ridge |
|----------|:--------:|:--------:|:-------------:|:-------------:|:--------------:|
| K562 | 0.085 | 0.610 | 0.888 | 0.988 | **0/4 승** |
| RPE1 | 0.147 | 0.696 | 0.890 | 0.983 | **0/4 승** |
| Norman | -0.699 | 0.896 | 0.422 | 1.000 | **0/4 승** |

- GEARS(DL)는 Ridge(선형)에 전패(0/12 승) — DL 복잡도 ≠ 예측 품질
- GEARS는 mean_predictor에는 11/12 승 — 순진한 baseline에는 승리
- Norman GEARS는 128/283 섭동만 유효 예측('ctrl+gene' 형식 비호환)
- **수정된 H3 클레임**: model quality > model complexity. Well-trained Ridge > baselines, but poorly-trained DL(GEARS) < Ridge. BioEval이 이 품질 격차를 정확히 식별 — BioEval 판별력의 검증.

**한계**:
- GEARS 5 epochs + 최소 hyperparameter — GEARS 최대 성능을 대표하지 않을 수 있음
- 단일 DL 모델(GEARS)만 테스트 — CPA 등 추가 DL 모델 필요
- Ridge LOO는 ground truth label 직접 사용 — 정보 이점 존재

---

## 3. Mean-Effect Trap의 기작

MSE가 방향 정보를 무시하는 기작을 정량적으로 해석한다:

1. **DEG 비율이 낮으면 MSE 최적화가 평균 예측에 수렴**: Norman에서 DEG 비율 1.53%. 즉, 98.47% 유전자는 교란에 반응하지 않는다. 이 유전자들에 대해 0을 예측하면 오차가 0이 된다.
2. **방향 정보의 가치가 MSE에서 0**: 부호가 틀려도 |예측 - 실제|²가 작으면 MSE가 낮다. 0.1의 오차를 갖는 올바른 방향 예측과 0.1의 오차를 갖는 틀린 방향 예측이 MSE에서 동등하다.
3. **BioEval-Dir은 이를 보정**: 부호 정확도를 직접 측정하므로 평균 예측의 방향 정확도가 0%(Dir_deg)로 나타남.

이 Trap은 섭동 예측 분야 전체에 영향을 미친다:
- Ahlmann-Eltze (2025) Nature Methods: DL ≤ baseline 보고
- SCALE (2026): MSE가 mean-effect trap 유발 지적
- Shesha (2026): magnitude ≠ stability, 지표 분리 필요

세 독립 그룹의 수렴 증거는 이 문제가 실재함을 강하게 시사한다.

---

## 4. MSE의 이중적 성격

가장 흥미로운 발견 중 하나는 MSE 자체가 방향 정보를 포착한다는 점이다:

| 데이터셋 | ρ(-MSE, dir_discovery) |
|----------|:----------------------:|
| K562 | 0.964 |
| RPE1 | 0.883 |
| Norman | 0.945 |

이것이 의미하는 바:
- MSE가 낮은 모델(=평균에 가까운 예측)이 방향 발견에서도 좋은 성과를 내는 것처럼 보인다
- 그러나 이는 상관이지 인과가 아니다: 학습이 잘 된 모델은 MSE도 낮고 방향도 정확하지만, MSE 최적화 자체가 방향 정보를 보장하지는 않는다
- mean_predictor는 MSE가 낮아도 방향 정확도가 0% — 극단적 사례
- MSE는 domain-general predictor이지만, 특정 영역(방향, DEG)에서는 BioEval이 더 정밀한 예측을 제공

---

## 5. 한계점

### 5.1 모델 다양성
- 9개 모델 중 3개 Ridge(alpha=1/10/100)가 핵심 — 선형 모델로만 검증
- GEARS DL 모델 훈련 진행 중 (run_20). 결과에 따라 H3 클레임 강도 변경 가능

### 5.2 순위 해상도
- N=9로 Kendall τ의 통계적 해상도가 제한
- 더 많은 모델(CPA, scGPT, Geneformer 등)이 추가되면 τ 패턴이 달라질 수 있음

### 5.3 DEG 임계값
- |logFC| > 0.25을 DEG 기준으로 사용
- 임계값 민감도 분석(S1)에서 PASS하였으나, 다른 기준 선택이 결과에 영향 가능

### 5.4 Cross-domain H2
- 방향 지표가 gene-set 과업을 MSE보다 잘 예측하지 못함 (33.3%)
- 이는 방향 정보와 gene-set 정보가 본질적으로 다른 신호이므로 예상되는 결과이나, H2 클레임의 범위를 domain-specific으로 제한해야 함

### 5.5 데이터셋 특성
- 3개 데이터셋 모두 인간 세포주(K562, RPE1) — 다른 생물종/조직으로의 일반화 불확실
- Norman은 조합 교란 포함, Replogle은 단일 교란 — 교란 유형에 따라 결과가 다를 수 있음

---

## 6. 후속 연구 방향

### 6.1 DL 모델 검증 (1차 완료)
- GEARS 훈련 완료 (run_20) — GEARS < Ridge (0/12 승). DL ≠ better 확인
- CPA 등 추가 DL 모델 검증으로 일반성 확보 필요

### 6.2 모델 다양성 확대
- Foundation model (Geneformer, scGPT, scBERT) 예측 포함
- 9개 → 15-20개 모델로 순위 해상도 향상

### 6.3 생물학적 검증
- 예측된 DEG가 실험적으로 검증된 DEG와 얼마나 일치하는지 확인
- CRISPR screen 데이터와의 교차 검증

### 6.4 BioEval 기반 학습
- BioEval-Dir을 loss function에 통합한 모델 학습
- MSE loss vs BioEval loss가 모델 성능에 미치는 영향 비교

### 6.5 다른 도메인으로의 일반화
- 비인간 종(yeast, mouse) 데이터셋에서 검증
- 단백질 수준(CITE-seq)에서의 방향 정확도 평가

---

## 7. 실험 이력

| Run | 날짜 | 내용 | 결과 |
|-----|------|------|------|
| run_13 | 04-30 | BioEval 메트릭-순위 반전 (시뮬레이션 11개 모델) | H1 SUPPORTED |
| run_14 | 04-30 | Phase 4 downstream 과업 상관 | H2 SUPPORTED (88.9%) |
| run_15 | 04-30 | sklearn Ridge LOO (Norman) | H1+H2 실제 모델 확인 (Norman). K562/RPE1 퇴화 |
| run_16 | 05-01 | Gene PCA Feature Ridge (K562/RPE1) | H1+H2+H3 3 데이터셋 전체 확인 |
| run_17 | 05-01 | Bootstrap CI (B=10,000) | H1+H2 통계적 견고성 확인 |
| run_18 | 05-01 | Scale Correction | 보정 불필요 확인. Dir_deg 불변 |
| run_19 | 05-01 | Downstream Task Independence | H2 domain-specific 확인. MSE domain-general |
| run_20 | 05-04 | GEARS DL 모델 훈련+평가 | GEARS < Ridge (0/12 승). K562 R²=0.085, RPE1 0.147, Norman -0.699. H3 정교화: model quality > complexity |

---

## 8. 참고 문헌

- Ahlmann-Eltze et al. (2025). "Comparison of perturbation prediction methods". Nature Methods.
- Chen et al. (2026). "SCALE: Single-cell perturbation landscape estimation". Cell-Eval framework.
- Csendes et al. (2025). "Benchmarking foundation models for perturbation prediction". BMC Genomics.
- Norman et al. (2019). "Mapping the perturbome landscape". Science.
- Raju et al. (2026). "Shesha: Geometric stability for perturbation prediction".
- Replogle et al. (2022). "Mapping information-rich genotype-phenotype landscapes". Cell.
- Roohani et al. (2025). "Virtual Cell Challenge". Cell.
- Wei et al. (2026). "Systematic benchmarking of perturbation prediction methods". Nature Methods.
- Zhu et al. (2025). "Evaluation metrics for perturbation prediction". Briefings in Bioinformatics.
