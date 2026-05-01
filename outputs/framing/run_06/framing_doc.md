# BioEval Framing — Run 06

**Date:** 2026-04-30 | **Direction:** Evaluation Metrics for Perturbation Prediction (Gap 2)
**Predecessor:** Literature Review run_07 (Gap Analysis) + Deep Dive

---

## 1. 연구 질문

**"섭동 예측 평가 지표가 생물학적 충실도를 측정하는가? 생물학적 충실도를 측정하는 지표를 설계하고, 이 지표 하에서 Ahlmann-Eltze의 'DL ≤ baseline' 위기가 해소되는가?"**

### 하위 질문

**RQ1 (지표-생물학 상관 진단):** 기존 평가 지표(MSE, R², Pearson, DEG overlap, PDS) 중 어떤 것이 생물학적 유용성(DEG 회복, 방향 정확도, downstream 과업 성과)과 상관하는가?

**RQ2 (지표 설계):** 생물학적으로 기반하고, 유전자 수준 분해능을 가지며, 방향을 인식하고, 구현에 견고한 평가 지표를 설계할 수 있는가? 이 지표가 기존 지표보다 downstream 생물학적 유용성을 더 잘 예측하는가?

**RQ3 (베이스라인 위기 해소):** 생물학적 충실도 지표 하에서 DL 모델이 simple baseline을 능가하는 체계(synergistic GI, 교세포 전이, 고다기능 섭동)가 존재하는가? 아니면 베이스라인 위기가 지표와 무관하게 실재하는가?

---

## 2. 신규성 사전 검사 (Novelty Safeguard)

### 2.1 경쟁 밀도 평가

| 용어 체계 | 직접 경쟁자 | 수준 |
|-----------|-----------|------|
| "perturbation prediction evaluation metric" | Zhu et al. (2025) AUPRC | 부분 — DEG 식별만, 방향/안정성/효과크기 없음 |
| "perturbation prediction benchmark" | Wei et al. (2026) Nature Methods 27방법×29데이터×6지표 | 부분 — 기존 지표로 벤치마크, 새 지표 설계 없음 |
| "foundation model perturbation benchmark" | Csendes et al. (2025) BMC Genomics | 부분 — FM vs baseline 비교만, 지표 분석 없음 |
| "Cell-Eval evaluation framework" | SCALE/Chen et al. (2026) PDCorr+DE overlap | 부분 — 다지표 프레임워크, 단 구현 민감성+이론적 분석 없음 |
| "geometric stability perturbation" | Shesha/Raju (2026) | 부분 — 단일 지표(안정성), 통합 프레임워크 없음 |

**총 직접 경쟁자: 5개 (부분적 겹침)**

⚠ 경쟁 밀도 ≥ 3 — 차별화 포인트 명확화 필요

### 2.2 차별화 분석

| 우리 기여 | 누가 이미 했는가? | 차별화 포인트 |
|-----------|-----------------|-------------|
| DEG 회복 평가 | Zhu (AUPRC), SCALE (DE overlap) | 우리: DEG precision+recall 곡선(단일 점 아님), 방향 정확도 결합 |
| 방향 일관성 | Shesha (안정성), SCALE (PDCorr) | 우리: 유전자 수준 방향 정확도(부호+크기), 안정성과 분리 |
| 효과크기 보정 | 없음 | 우리: 예측 logFC vs 실제 logFC의 보정 곡선 — 과소/과대 예측 체계적 탐지 |
| 베이스라인 위기 진단 | Ahlmann-Eltze (DL≤baseline), Csendes (FM≤mean) | 우리: **지표 교체 시 순위 반전 여부의 정량적 분석** — 이것이 핵심 차별화 |
| 통합 프레임워크 | SCALE (Cell-Eval) | 우리: 이론적 속성 분석(지표가 무엇을 측정하는가?), 구현 견고성 테스트 |
| 지표-downstream 과업 상관 | 없음 | 우리: 지표 순위 → AL/Hit prioritization 성과 예측력 — 완전 신규 |

### 2.3 기여 분해 테스트

| 기여 성분 | 독립 신규? | 비고 |
|-----------|----------|------|
| DEG precision/recall 곡선 | 부분 | AUPRC(Zhu)와 겹침. 차이: 곡선 전체 vs 단일 수치, 방향 결합 |
| 유전자 수준 방향 정확도 | 신규 | Shesha는 섭동 수준 안정성, 우리는 유전자×섭동 수준 부호 정확도 |
| 효과크기 보정 분석 | 신규 | logFC 예측 vs 실제의 체계적 편향 탐지 — 누구도 안 함 |
| 지표-순위 반전 정량화 | 신규 | MSE 순위 vs BioEval 순위의 Kendall τ, 순위 변동 유전자/섭동 특성 분석 |
| 지표-downstream 과업 상관 | 신규 | 지표 선택이 AL 효율, hit prioritization 정밀도에 미치는 영향 — 완전 미탐색 |
| 최소 보고 기준 제안 | 부분 | DOME(Walsh 2020)이 일반 ML, 우리는 섭동 예측 특화 |

**독립 신규 성분: 3/6 (50%)** — 조합 신규성 포함하면 4/6. 위험하지만 치명적이지 않음.

### 2.4 실패 모드 분류

| 실패 모드 | 징후 있음? | 근거 |
|-----------|----------|------|
| 포화 시장 | ⚠ 주의 | 5개 부분 경쟁자. 단 직접 동일 접근(지표 설계+위기 해소+downstream 상관)은 0개 |
| 사소한 개선 | 낮음 | RQ3(순위 반전)은 이진 결과 — 반전되거나 안 되거나. 둘 다 중요 |
| 재포장 | 낮음 | Wei et al.은 벤치마크(기존 지표 사용), 우리는 지표 설계(새 지표 제안). 본질적으로 다른 과업 |
| 평가 함정 | ⚠ 주의 | BioEval 자체의 구현 민감성 위험. 이것을 메타-분석으로 방어해야 |
| 데이터 의존 | 낮음 | Ahlmann-Eltze 7+ 벤치마크 + Replogle + Norman으로 교차 검증 |

**2개 실패 모드에서 ⚠ — 방어 전략 필요**

---

## 3. 핵심 차별화 전략: "Metric-Ranking Reversal Analysis"

Wei et al. (2026)과 Zhu et al. (2025)가 이미 존재하므로, **순수 벤치마크나 단일 지표 제안으로는 차별화 부족**. 우리의 핵심 차별화는:

### 3.1 핵심 질문: "지표를 바꾸면 모델 순위가 바뀌는가?"

이것이 Wei et al.이나 Zhu et al.이 묻지 않은 질문이다:
- Wei et al.: 6개 지표로 벤치마크 → "어떤 모델이 좋은가?"
- Zhu et al.: AUPRC 제안 → "R² 대신 AUPRC를 쓰라"
- 우리: **"지표 선택에 따라 모델 순위가 어떻게 변하는가? 어떤 지표가 downstream 생물학적 유용성을 예측하는가?"**

이것은 메타-평가(meta-evaluation) — 지표 자체를 평가하는 것이다.

### 3.2 구체적 분석

1. **Kendall τ 분석**: MSE 순위 vs BioEval 순위의 순위 상관. τ < 0.5이면 지표가 순위를 바꿈
2. **순위 변동 분석**: 어떤 모델/섭동/유전자에서 순위가 바뀌는가? (예: synergistic GI에서만 DL이 우위?)
3. **지표-downstream 상관**: 각 지표 순위 vs AL 효율 순위, hit prioritization 순위. 어느 지표가 downstream을 예측하는가?
4. **효과크기 보정**: 예측 logFC의 체계적 과소/과대 예측 탐지. MSE는 이것을 놓침

---

## 4. 평가 전략

### 4.1 RQ1: 지표-생물학 상관 진단

| 지표 | 범주 | 측정 대상 |
|------|------|----------|
| MSE | 분포 매칭 | 예측-관측 평균제곱오차 |
| R² | 분포 매칭 | Pearson 상관의 제곱 |
| Pearson | 분포 매칭 | 선형 상관 |
| DE overlap | DEG 회복 | DEG 집합 교집합/합집합 |
| AUPRC (Zhu) | DEG 회복 | DEG precision-recall 곡선 아래 면적 |
| PDCorr (SCALE) | 방향 | 섭동 방향 상관 |
| Shesha stability | 방향 | 기하학적 안정성 |
| **BioEval-Dir** (제안) | 방향 | 유전자×섭동 수준 부호 정확도 + 크기 비율 |
| **BioEval-Cal** (제안) | 보정 | logFC 예측 vs 실제 보정 곡선 |
| **BioEval-Composite** (제안) | 통합 | Dir + Cal + DEG 가중합 |

**downstream 과업 (지표가 예측해야 할 것):**
1. AL 효율: 지표로 선택된 모델이 AL에서 더 나은가?
2. Hit prioritization: 지표로 선택된 모델이 top-k DEG를 더 잘 찾는가?
3. 교세포 전이: 지표로 선택된 모델이 다른 세포유형에서 더 잘 작동하는가?

### 4.2 RQ2: BioEval 지표 설계

**BioEval-Dir (유전자 수준 방향 정확도):**
- 각 유전자 g, 섭동 p에 대해: sign(ŷ_gp) == sign(y_gp)이면 +1, 아니면 0
- 크기 비율: |ŷ_gp| / |y_gp| (과소 < 1, 과대 > 1)
- 집계: 방향 정확도 = 평균(부호 일치), 크기 비율 = 중앙값(|ŷ/y|)

**BioEval-Cal (효과크기 보정):**
- logFC 예측 vs 실제의 산점도 → 회귀 기울기 (이상=1.0)
- 과소 예측 편향: 기울기 < 1.0인 유전자/섭동 비율
- 과대 예측 편향: 기울기 > 1.0인 유전자/섭동 비율

**BioEval-Composite:**
- 방향 정확도 × DEG precision × 보정 기울기의 가중 조합
- 가중치는 downstream 과업 상관에서 학습

### 4.3 RQ3: 베이스라인 위기 해소

**실험 설계:**
1. Ahlmann-Eltze 벤치마크의 모델 예측을 확보 (또는 재현)
2. MSE 순위 계산 → Baseline 순위 확인
3. BioEval 지표로 재평가 → 순위 변화 측정
4. 순위 변동이 체계적인지 확인:
   - synergistic GI에서만 DL이 우위?
   - 교세포 전이 체계에서만 DL이 우위?
   - 고다기능 섭동에서만 DL이 우위?

**성공 기준:**
- 순위 반전(Kendall τ < 0.5): BioEval이 MSE와 다른 모델을 선택 → 지표가 원인이었음
- 순위 유지(Kendall τ > 0.7): BioEval 하에서도 baseline이 우위 → 위기 실재
- 부분 반전(0.5 ≤ τ ≤ 0.7): 특정 체계에서만 DL 우위 → 미묘한 결과

**둘 다 발표 가능:**
- 순위 반전 → "베이스라인 위기는 지표 아티팩트" (높은 임팩트)
- 순위 유지 → "베이스라인 위기는 실재, BioEval은 더 나은 진단 도구" (중간 임팩트)
- 부분 반전 → "DL은 특정 체계에서만 우위, BioEval이 이를 드러냄" (높은 임팩트)

---

## 5. 베이스라인 계층

1. **우연 수준**: 랜덤 예측 (MSE 최대, BioEval-Dir=0.5, AUPRC=π₀)
2. **Mean predictor**: 관측치 평균 (Ahlmann-Eltze 베이스라인)
3. **Additive linear**: Y = GW^T P + b (Ahlmann-Eltze 최우수)
4. **scGPT + linear**: 선형 모델 + 사전학습 임베딩 (Ahlmann-Eltze에서 최고)
5. **CPA**: 조합 오토인코더 (표준 비교)
6. **GEARS**: GNN + GRN (조합 예측 표준)
7. **SCALE**: LLaMA 기반 파운데이션 모델 (최신 FM)

---

## 6. 데이터

| 데이터셋 | 용도 | 비고 |
|----------|------|------|
| Replogle 2022 | RQ1-3 주 평가 | K562+RPE1, 848 공유 섭동, 교세포 |
| Norman 2019 | RQ1-3 조합 평가 | 128 double-KO, GI ground-truth |
| Ahlmann-Eltze 벤치마크 | RQ3 직접 재현 | 7+ 벤치마크, 모델 예측 확보 필요 |
| PBMC (Zhu 2025) | RQ1 AUPRC 비교 | 7 cell types, IFN-γ 자극 |
| PORTAL 2026 | RQ3 대규모 검증 | 665K pairwise (선택) |

---

## 7. 타겟 성능 수치

| RQ | 지표 | 베이스라인 | 타겟 | 근거 |
|----|------|----------|------|------|
| RQ1 | 지표-AL 상관 (Spearman) | 0 (MSE-AL 상관) | > 0.5 | MSE가 AL을 예측 못하면 BioEval이 해야 |
| RQ1 | 지표-DEG 상관 (Spearman) | AUPRC 기준치 | BioEval > MSE by ≥0.1 | Zhu: R²-AUPRC 불일치 이미 입증 |
| RQ2 | BioEval-Dir 방향 정확도 | 0.5 (우연) | > 0.7 | 생물학적 방향 일치가 70%+이어야 유용 |
| RQ2 | BioEval-Cal 보정 기울기 | — | 0.8-1.2 범위 | 과소/과대 예측 20% 이내 |
| RQ3 | Kendall τ (MSE vs BioEval 순위) | 1.0 (동일) | < 0.5 (반전) 또는 > 0.7 (유지) | 이진 결과 — 둘 다 의미 있음 |
| RQ3 | Synergistic GI에서 DL>baseline 비율 | 0% (MSE 하) | > 30% (BioEval 하) | 특정 체계에서 DL 우위 탐지 |

---

## 8. 위험과 대응

| 위험 | 확률 | 영향 | 대응 |
|------|------|------|------|
| Wei et al. (2026)과 중복 | 중간 | 높음 | 우리는 지표 설계+순위 반전 분석, Wei는 기존 지표 벤치마크. 본질적 차이 |
| BioEval도 구현 민감성 | 중간 | 중간 | SCALE이 지적한 문제. 임계값 스윕으로 견고성 테스트 |
| 순위 반전이 안 일어남 | 중간 | 낮음 | 부정 결과도 발표 가능. "위기 실재" = 중요 발견 |
| Ahlmann-Eltze 예측 확보 불가 | 낮음 | 높음 | 코드 공개 확인 필요. 미공개 시 Replogle+Norman으로 재현 |
| 생물학적 유용성의 정의 모호 | 중간 | 중간 | DEG 회복, AL 효율, hit prioritization 3가지로 구체화 |

---

## 9. 기존 FCR-ICM 프로젝트에서의 교훈 적용

FCR-ICM의 실패 원인이 이 프로젝트에서 재발하지 않도록:

1. **경쟁자 조기 확인**: FCR-ICM은 8 runs 후 BuDDI/C3TL 발견. BioEval은 framing 단계에서 5개 부분 경쟁자 이미 확인
2. **ICM 부정 결과**: ICM이 에피스태시스 전이에 기여 못함(0.966×). 이것은 "단순 모델이 복잡한 것보다 낫다"는 같은 맥락 — 베이스라인 위기의 국소적 증거
3. **CPA>FCR (0.430 vs 0.367)**: 이것도 "단순 모델이 낫다"의 증거. 하지만 MSE 기준이었음 — BioEval에서는 반전 가능?
4. **지표 선택이 결론 변경**: prod rho=0.437 PASS vs A7 rho=0.326 PARTIAL — 같은 데이터에서 지표가 결론을 바꿈. 이것이 RQ3의 직접적 동기

---

## 10. 핵심 논문

| 논문 | DOI/arXiv | 역할 | 위협 수준 |
|------|-----------|------|----------|
| Ahlmann-Eltze et al. (2025) | 10.1038/s41592-025-02772-6 | 흡연총: DL ≤ baseline | 낮음 (우리의 출발점) |
| Wei et al. (2026) | 10.1038/s41592-025-02980-0 | 27방법 벤치마크 | **높음** — 직접 경쟁 |
| Zhu et al. (2025) | 10.1093/bib/bbaf426 | AUPRC 지표 | 중간 — 단일 지표만 |
| Csendes et al. (2025) | 10.1186/s12864-025-11600-2 | FM ≤ mean | 낮음 — 베이스라인 위기 확인 |
| SCALE/Chen et al. (2026) | arXiv:2603.17380 | PDCorr+Cell-Eval | 중간 — 부분 겹침 |
| Shesha/Raju (2026) | arXiv:2604.16642 | 기하학적 안정성 | 낮음 — 단일 지표 |
| Roohani et al. (2025) | 10.1016/j.cell.2025.06.008 | Virtual Cell Challenge | 낮음 — PDS 벤치마크 |
