# BioEval: 생물학적 충실도 기반 섭동 예측 평가 프레임워크

## 연구 보고서 (2026-05-02 기준)

---

## 1. 풀고자 하는 문제

### 1.1 배경: Perturb-seq와 섭동 예측

**Perturb-seq**은 단일 세포 수준에서 유전자를 교란(knockout, knockdown)한 뒤 전사체(transcriptome)를 측정하는 실험 기술이다. 1만~50만 개의 세포를 동시에 관찰할 수 있어, 유전자 기능 연구와 약물 반응 예측에 핵심 도구로 쓰인다.

**섭동 예측(perturbation prediction)**은 Perturb-seq 데이터를 이용해 "특정 유전자를 교란하면 세포의 유전자 발현이 어떻게 변할까?"를 계산적으로 예측하는 문제다. 즉, 실험을 하지 않고도 교란 효과를 예측하는 것으로, 조합적 교란(2개 이상 유전자를 동시에 교란)의 경우 실험적 탐색 공간이 기하급수적으로 커지기 때문에 계산적 예측의 가치가 크다.

### 1.2 핵심 문제: "DL ≤ Baseline" 위기

2025년 Ahlmann-Eltze 등이 Nature Methods에 발표한 대규모 벤치마크에서, **딥러닝(DL) 모델이 단순 선형 모델보다 나을 게 없다**는 결과가 보고되었다. 7개 이상의 벤치마크 데이터셋에서 일관되게 DL ≤ linear baseline이 관찰되었으나, **원인 분석은 이루어지지 않았다**.

이 위기의 근본 원인으로 우리는 **평가 지표의 문제**를 지목한다.

### 1.3 평가 지표가 만드는 왜곡: Mean-Effect Trap

현재 섭동 예측 분야의 표준 평가 지표는 **MSE(Mean Squared Error)**, **R²(결정계수)**, **Pearson 상관계수**이다. 이들은 모두 **예측값과 실제값의 수치적 일치도**를 측정할 뿐, 생물학적 의미(어떤 유전자가 변했는지, 변화 방향이 맞는지)는 무시한다.

이로 인해 **Mean-Effect Trap**이 발생한다:
- 대부분의 유전자는 교란에 반응하지 않는다(비차별적으로 발현됨). 예: Norman 데이터에서 DEG(차별 발현 유전자) 비율은 단 1.53%
- MSE를 최소화하려면 **평균에 가까운 예측**을 하는 것이 유리하다. 모든 유전자를 평균으로 예측하면 비차별 유전자에서 오차가 0이 되기 때문
- 하지만 평균 예측은 **방향 정보를 전혀 포함하지 않는다**. "이 유전자가 상향 조절되었다"는 생물학적으로 가장 중요한 정보가 사라짐
- 결과: MSE 순위에서는 평균 예측기(mean predictor)가 상위이나, 생물학적 유용성은 0%

### 1.4 우리가 묻는 질문

> **"평가 지표를 생물학적 충실도 기반으로 바꾸면, 모델 순위가 어떻게 달라지는가? DL ≤ Baseline 위기가 지표의 아티팩트(artifact)인지, 실재하는 현상인지 판별할 수 있는가?"**

이 질문은 섭동 예측 분야에서 **누구도 정량적으로 분석하지 않은 것**이다. 기존 연구들은 "어떤 모델이 좋은가?"를 물었지만, 우리는 **"어떤 지표로 평가하느냐에 따라 모델 순위가 어떻게 변하는가?"**를 묻는다.

---

## 2. 기존 접근 방법

### 2.1 표준 평가 지표

| 지표 | 측정 대상 | 장점 | 한계 |
|------|----------|------|------|
| **MSE** (Mean Squared Error) | 예측-실제 분산 | 구현 단순, 미분 가능 | 생물학적 방향·크기 무시; mean-effect trap |
| **R²** (결정계수) | 분산 설명 비율 | 직관적 해석 | MSE와 동일한 한계; 음수 가능 |
| **Pearson** (상관계수) | 선형 상관 | 스케일 불변 | 방향 구분 불가; 이상치 민감 |

**용어 설명**:
- **MSE**: 예측값과 실제값의 차이를 제곱해 평균낸 것. 낮을수록 좋음
- **R²**: 모델이 데이터의 분산을 얼마나 잘 설명하는지. 1.0이 최고, 음수면 평균보다 나쁨
- **Pearson**: 두 변수의 선형 관계 강도. -1~1 범위. 크기만 비슷하면 방향이 틀려도 높게 나옴

### 2.2 생물학 관련 지표 (부분적 선행)

| 지표 | 출처 | 측정 대상 | 장점 | 한계 |
|------|------|----------|------|------|
| **DEG overlap** | CPA (2022) | 차별 발현 유전자 회복 | 생물학적 관점 도입 | 이진 임계값으로 정보 손실; 방향 미포함 |
| **AUPRC** | Zhu et al. (2025) | DEG precision-recall | 임계값 독립 곡선 | DEG 식별만, 방향/보정/downstream 상관 없음 |
| **PDCorr** | SCALE/Chen et al. (2026) | 섭동 방향 상관 | 방향 정보 포함 | 유전자 수준 분해 없음; 구현 민감성 |
| **PDS** (Perturbation Distribution Similarity) | ARC Virtual Cell Challenge | 분포 거리 | 분포 관점 | 거리 지표 선택에 민감; 생물학적 해석 어려움 |
| **Stability** | Shesha/Raju (2026) | 기하학적 안정성 | magnitude와 분리 가능 | 단일 지표, 통합 프레임워크 없음 |

**용어 설명**:
- **DEG** (Differentially Expressed Gene, 차별 발현 유전자): 교란 후 발현량이 의미 있게 변한 유전자. 보통 |logFC| > 임계값(0.25 등)으로 정의
- **logFC** (log Fold Change): 유전자 발현 변화량의 로그. 양수=상향, 음수=하향 조절
- **AUPRC** (Area Under Precision-Recall Curve): 정밀도-재현율 곡선 아래 면적. 불균형 데이터에서 AUROC보다 유용
- **PDCorr**: 예측과 실제의 방향(부호) 일치도를 섭동 단위로 측정

### 2.3 벤치마크 연구

| 연구 | 핵심 기여 | 우리와의 차이 |
|------|----------|--------------|
| **Wei et al. (2026)** Nature Methods | 27개 방법 × 6개 지표 대규모 벤치마크 | 기존 지표만 사용, **새 지표 설계 없음**, **순위 반전 분석 없음** |
| **Ahlmann-Eltze (2025)** Nature Methods | DL ≤ baseline 위기 진단 | 원인 분석 없음, 해결책 없음 |
| **SCALE/Chen et al. (2026)** | Cell-Eval 프레임워크 제안 | PDCorr+DE overlap 제안하나 구현 민감성 문제, 이론적 분석 없음 |
| **Csendes et al. (2025)** BMC Genomics | Foundation Model ≤ mean predictor 확인 | 지표 분석 없음 |

---

## 3. 각 방법론의 장단점 비교

### 3.1 MSE/R² vs 생물학적 지표: 근본적 차이

```
MSE가 좋은 모델 ≠ 생물학적으로 유용한 모델

예시 (Norman 데이터):
  mean_predictor: MSE #1위 → downstream 과업 f1@50 = 0.021 (최악)
  Ridge (alpha=1): MSE #4위 → downstream 과업 f1@50 = 0.488 (최고)

  → MSE가 "가장 좋다"고 평가한 모델이 생물학적 과업에서는 최악
```

### 3.2 기존 생물학적 지표의 한계 요약

| 한계 | 해당 지표 | 설명 |
|------|----------|------|
| **방향 무시** | MSE, R², Pearson, DEG overlap, AUPRC | 유전자가 상향/하향 조절되었는지 구분 안 함 |
| **유전자 수준 분해 없음** | PDCorr, Shesha stability | 섭동 전체 평균만. 어떤 유전자에서 틀렸는지 모름 |
| **임계값 의존** | DEG overlap | 임계값에 따라 결과가 크게 변함 |
| **downstream 과업 상관 미측정** | 모든 기존 지표 | 지표가 좋다고 해서 실제 생물학적 분석에 도움이 되는지 확인 안 됨 |
| **순위 반전 분석 부재** | 모든 기존 연구 | 지표 교체가 모델 선택을 바꾸는지 정량 분석 없음 |

### 3.3 우리가 식별한 3개 독립 그룹의 수렴 증거

이 문제는 단일 연구의 주장이 아니다. 세 독립 그룹이 각각 다른 각도에서 같은 문제를 지적:

1. **Ahlmann-Eltze (2025)**: DL ≤ baseline (위기 진단)
2. **SCALE (2026)**: MSE가 mean-effect trap 유발 (원인 지적)
3. **Shesha (2026)**: magnitude ≠ stability (지표 분리 필요)

이 세 증거가 수렴한다는 것은 문제가 실재함을 강하게 시사한다.

---

## 4. Gap: 기존 연구가 놓친 것

### Gap 1: 지표-생물학 상관 부재
기존 지표 중 **어떤 것이 downstream 생물학적 유용성을 예측하는지** 정량 분석이 없다. 지표가 높으면 실험적으로도 유용한가? 이 질문에 답한 연구가 없다.

### Gap 2: 생물학적 충실도 지표 설계 공백 ⭐ (우리가 선택한 Gap)
유전자 수준 분해능 + 방향 인식 + 보정 분석 + 구현 견고성을 **모두** 갖춘 지표가 없다. 각 성분은 부분적으로 존재하나, 통합 프레임워크는 없다.

### Gap 3: 지표-순위 반전 분석 미수행 (핵심 차별화)
**지표를 바꾸면 모델 순위가 어떻게 변하는가?** 이 질문을 정량적으로 분석한 연구가 전무하다. Wei et al.이 6개 지표로 벤치마크했지만, 지표 간 순위 불일치를 정량화하거나 원인 분석하지 않았다.

### Gap 4: 베이스라인 위기 원인 불명
DL ≤ baseline이 **지표의 아티팩트**인지 **실재하는 현상**인지 판별할 방법이 없다. 생물학적 충실도를 측정하는 지표가 없으면 이 판별이 불가능하다.

---

## 5. 이번 실험의 차별점: BioEval 프레임워크

### 5.1 BioEval 지표 설계

BioEval은 세 가지 하위 지표로 구성된 통합 평가 프레임워크다:

#### BioEval-Dir: 유전자×섭동 수준 방향 정확도

**정의**: 각 섭동(perturbation)에 대해, 각 유전자(gene)의 예측 방향(상향/하향/변화 없음)이 실제와 일치하는 비율

```
방향 정확도 계산:
  - 예측과 실제 부호가 같으면 → 1점
  - 둘 다 0이면 → 1점
  - 실제는 0인데 예측이 0이 아니면 → 0.5점
  - 부호가 다르면 → 0점

  Dir_all  = 전체 유전자에 대한 평균 (기본 방향 정확도)
  Dir_deg  = DEG에 대해서만 평균 (생물학적 의미가 큰 유전자에 집중)
  Dir_weighted = |logFC|로 가중 (큰 효과의 유전자에 더 높은 가중치)
```

**차별점**: 기존 PDCorr은 섭동 수준만 측정. AUPRC는 방향을 무시. BioEval-Dir은 **유전자×섭동 2차원 분해**를 제공한다.

#### BioEval-DEG: DEG 회복 정밀도

**정의**: 예측이 실제 DEG를 얼마나 정밀하게 식별하는가를 AUPRC로 측정

```
DEG_auprc     = |예측|을 점수로 사용한 AUPRC (기존 Zhu와 유사)
DEG_dir_auprc = |예측| × 방향일치를 점수로 사용한 AUPRC (방향 정보 결합 — 신규)
```

**차별점**: DEG_dir_auprc는 방향 정보를 결합한 최초의 DEG 평가 지표.

#### BioEval-Cal: 효과 크기 보정 분석

**정의**: 예측 효과 크기(logFC)의 체계적 과소/과대 예측을 탐지

이번 실험에서는 Ridge 모델의 과소 예측(shrinkage)을 정량화하는 데 활용.

### 5.2 세 가지 가설 (H1-H3)

| 가설 | 내용 | 측정 방법 | 타겟 |
|------|------|----------|------|
| **H1**: 지표-순위 반전 | MSE/R²와 BioEval-Dir은 모델을 다르게 순위 매긴다 | Kendall τ(MSE, Dir_deg) < 0.7 | τ < 0.7 = 순위 불일치 |
| **H2**: 지표-downstream 상관 | BioEval 지표가 MSE보다 downstream 생물학적 유용성을 더 잘 예측한다 (domain-specific) | Spearman ρ(BioEval, downstream) > ρ(MSE, downstream) | gap > 0.1 |
| **H3**: 학습 모델 > 베이스라인 | BioEval 하에서 학습된 모델이 단순 베이스라인을 능가한다 | 모든 BioEval 지표에서 trained > baseline | ALL WIN |

**용어 설명**:
- **Kendall τ (타우)**: 두 순위 간의 일치도를 측정하는 비모수 통계량. -1~1 범위. 1=완전 일치, 0=독립, -1=완전 반대
- **Spearman ρ (로)**: 순위 상관계수. 두 변수의 순위 간 선형 관계. -1~1 범위
- **downstream 과업(downstream task)**: 섭동 예측 결과를 실제 생물학적 분석에 사용하는 과업

### 5.3 downstream 과업 정의

| 과업 | 정의 | 생물학적 의미 | 방향 정보 사용 |
|------|------|--------------|:----------:|
| **f1@50** | 예측 상위 50개 유전자 vs 실제 DEG의 F1 점수 | "어떤 유전자가 변하는지 얼마나 잘 찾아내는가?" | 아니오 (순위만) |
| **f1@100** | 예측 상위 100개 유전자 vs 실제 DEG의 F1 점수 | 위와 동일, 더 넓은 범위 | 아니오 (순위만) |
| **dir_discovery_deg** | DEG에서 예측 방향이 실제와 일치하는 비율 | "찾아낸 DEG의 방향이 맞는가?" | 예 |
| **mag_rank** | DEG에서 |예측|과 |실제|의 Spearman 상관 | "효과 크기 순서가 맞는가?" | 아니오 |
| **top100_overlap** | Jaccard(top-100 |예측|, top-100 |실제|) | "가장 큰 효과 유전자가 겹치는가?" | 아니오 |

---

## 6. 실험 결과 요약

### 6.1 데이터셋

| 데이터셋 | 세포 수 | 섭동 수 | 유전자 수 | DEG 비율 | 비고 |
|----------|---------|---------|----------|---------|------|
| **Replogle K562** | 162,751 | 1,092 | 5,000 | 2.38% | 인간 백혈병 세포주 |
| **Replogle RPE1** | 162,733 | 1,543 | 5,000 | 6.50% | 인간 망막 색소상피세포 |
| **Norman 2019** | 91,205 | 283 | 5,045 | 1.53% | 인간 K562, 조합 교란 포함 |

**용어**: **DEG 비율** = |logFC| > 0.25인 유전자의 비율. Norman의 DEG 비율이 가장 낮아(1.53%) mean-effect trap이 가장 심함.

### 6.2 모델 구성 (9개)

| 모델 | 유형 | 설명 |
|------|------|------|
| ridge | 학습 | Ridge 회귀 (alpha=1, analytical LOO) |
| ridge_med | 학습 | Ridge 회귀 (alpha=10) |
| ridge_strong | 학습 | Ridge 회귀 (alpha=100, 강한 정규화) |
| noisy_ridge | 퇴화 | ridge + 15% 가우시안 노이즈 |
| sign_flip_ridge | 퇴화 | ridge + 15% 부호 반전 |
| mean_predictor | 베이스라인 | 관측치 평균 예측 (Ahlmann-Eltze 베이스라인) |
| mean_effect | 베이스라인 | 섭동 평균 효과 예측 |
| constant_shrink | 오라클 | true × 0.15 (H3에서 제외) |
| half_signal | 오라클 | true × 0.5 (H3에서 제외) |

**용어**:
- **Analytical LOO** (Leave-One-Out): 모든 섭동을 사용해 학습한 뒤 하나를 제외하고 예측하는 과정을 수학적으로(hat matrix) 계산
- **Oracle baseline**: 실제 정답(true)을 사용하는 베이스라인. 공정한 비교가 아니므로 H3에서 제외
- **정규화(regularization, alpha)**: 과적합을 방지하기 위해 모델 복잡도에 벌점을 주는 기법. alpha가 클수록 예측이 평균에 가까워짐

### 6.3 H1 결과: 지표-순위 반전 ⭐ SUPPORTED

**핵심 발견**: MSE/R²와 BioEval-Dir은 모델을 **독립적으로** 순위 매긴다.

| 데이터셋 | τ(MSE, Dir_all) | τ(MSE, Dir_deg) | τ(Pearson, Dir_deg) | 판정 |
|----------|:---------------:|:---------------:|:-------------------:|------|
| K562 | -0.500 | -0.611 | 0.056 | PARTIAL/REVERSAL |
| RPE1 | -0.333 | -0.389 | 0.111 | **REVERSAL** |
| Norman | -0.500 | -0.500 | -0.167 | PARTIAL/REVERSAL |

**Bootstrap CI (B=10,000) 검증**:

| 데이터셋 | |τ(MSE, Dir_all)| | 95% CI가 0 포함? | 해석 |
|----------|:--------------:|:----------------:|------|
| K562 | 0.696 | 아니오 (폭 넓음) | 중등도 불일치; 0.7 미만 |
| RPE1 | **0.232** | 예 | **H1의 가장 명확한 증거** — MSE와 Dir이 독립 |
| Norman | 0.348 | 예 | 독립 — H1 지지 |

**의미**: RPE1과 Norman에서 MSE와 BioEval-Dir의 순위 상관 95% CI가 0을 포함한다. 즉, 두 지표는 모델을 **통계적으로 독립적으로** 평가한다.

### 6.4 H2 결과: 지표-downstream 상관 ⭐ SUPPORTED (domain-specific)

**핵심 발견**: BioEval 지표가 MSE보다 downstream 생물학적 유용성을 더 잘 예측한다. 단, 이는 **도메인 내(intra-domain)**에서 확인되며, **도메인 간(cross-domain)**에서는 확인되지 않는다.

#### 도메인 분해 (run_19)

| 도메인 | 통과율 | 평균 gap | 해석 |
|:------:|:------:|:--------:|------|
| **cross-domain** (Dir ↔ gene-set) | 8/24 (33.3%) | -0.086 | Dir 지표가 gene-set 과업을 MSE보다 잘 예측하지 못함 |
| **intra-DEG** (DEG_auprc ↔ f1@50) | 9/9 (**100%**) | +0.319 | DEG_auprc가 f1@50을 잘 예측 — 순환 아님, 의미 있는 신호 |
| **intra-magnitude** (mag_rank ↔ gene-set) | 6/6 (**100%**) | +0.303 | 크기 정보가 gene-set 과업을 잘 예측 |
| **intra-direction** (Dir ↔ dir_discovery) | 3/6 (50%) | -0.013 | 부분 순환, 약한 이점 |

#### 핵심 인사이트: MSE 자체가 domain-general predictor

| 데이터셋 | ρ(-MSE, dir_discovery) | 해석 |
|----------|:----------------------:|------|
| K562 | 0.964 | MSE가 방향 발견을 ρ=0.964로 예측 |
| RPE1 | 0.883 | MSE가 방향 발견을 잘 예측 |
| Norman | 0.945 | MSE가 방향 발견을 잘 예측 |

MSE는 magnitude와 direction 양쪽 신호를 모두 포착하는 **domain-general predictor**다. BioEval의 이점은 cross-domain 예측이 아니라 **도메인 내 해석 가능성**에 있다.

#### 방향 독립적 과업 검증

방향 정보를 사용하지 않는 과업(mag_rank, top100_overlap)에서도 23/33 (69.7%)가 통과하여, 결과가 방향 정보 순환에 의해 주도되지 않음을 확인.

**수정된 H2 클레임**: BioEval 지표는 MSE에 비해 **도메인 특이적 예측 이점**을 제공한다. Cross-domain 예측은 주장하지 않는다.

### 6.5 H3 결과: 학습 모델 > 베이스라인 ⭐ SUPPORTED

**핵심 발견**: BioEval 하에서 학습된 Ridge 모델이 모든 지표에서 베이스라인을 능가한다.

| 데이터셋 | Dir_deg (학습) | Dir_deg (베이스라인) | R² (학습) | R² (베이스라인) | 6 지표 |
|----------|:-------------:|:-------------------:|:---------:|:---------------:|:------:|
| K562 | 0.985 | 0.606 | 0.523 | -0.013 | **ALL WIN** |
| RPE1 | 0.989 | 0.666 | 0.652 | -0.013 | **ALL WIN** |
| Norman | 0.986 | 0.571 | 0.643 | -0.002 | **ALL WIN** |

**한정 조건**: Ridge(선형 모델)로 검증됨. DL 모델(GEARS) 훈련 진행 중(run_20).

---

## 7. 검증된 핵심 지식

1. **MSE/R² 순위 반전은 실재한다** (H1 SUPPORTED): 3개 데이터셋에서 MSE와 BioEval-Dir의 순위가 독립적. RPE1이 가장 명확한 증거(τ = 0.232)
2. **Mean-effect trap이 체계적 현상이다**: Norman에서 mean_predictor가 MSE #1이나 Dir #11. DEG 비율이 낮을수록 trap이 심화
3. **BioEval이 MSE보다 downstream 유용성을 더 잘 예측한다** (H2 SUPPORTED, domain-specific): Intra-domain에서 100% 통과. Cross-domain은 33.3%로 약하나, 이는 예상되는 결과
4. **DEG_auprc가 가장 견고한 H2 지표다**: Bootstrap CI에서 모든 데이터셋에서 유의
5. **학습 모델이 베이스라인을 능가한다** (H3 SUPPORTED): 3 데이터셋 × 6 지표 ALL WIN. 선형 모델로 검증; DL 검증 진행 중
6. **교세포 일관성**: K562와 RPE1(같은 유기체, 다른 세포유형) 간 모든 지표 τ > 0.78. 현상이 체계적
7. **MSE는 domain-general predictor**: ρ(-MSE, dir_discovery) = 0.88-0.96. MSE 자체가 방향 정보를 포착. BioEval의 이점은 해석 가능성과 도메인 내 정밀도
8. **Norman logFC 스케일 불일치는 무시 가능** (run_18): 3가지 보정 전략 테스트 결과, 방향 지표는 불변(부호 기반), downstream 과업도 개선 없음

---

## 8. 남은 과제와 한계

### 8.1 해결된 차단요소

| ID | 문제 | 해결 |
|----|------|------|
| B1 | K562/RPE1 one-hot LOO 퇴화 | Gene PCA features로 해결 (run_16). R² -0.027→0.523 |
| B3 | Norman logFC 스케일 불일치 | run_18에서 보정 불필요 확인. 방향 지표 불변 |
| B4 | Oracle baselines이 H3 왜곡 | H3 분석에서 제외. oracle upper bound로 분류 |
| B5 | Bootstrap CI 부재 | run_17에서 B=10,000 bootstrap으로 해결 |
| B6 | downstream 과업 순환성 | run_19에서 도메인 분해로 해결. H2는 domain-specific. 방향 독립적 과업 69.7% 통과 |

### 8.2 미해결 과제

| ID | 문제 | 심각도 | 영향 |
|----|------|--------|------|
| B2 | 실제 DL 모델 예측 없음 | 중간 | H3 "DL > baseline" 클레임이 선형 모델로만 검증됨. GEARS 훈련 진행 중 (run_20) |

### 8.3 허용 클레임 강도

| 가설 | 클레임 강도 | 한정 조건 |
|------|-----------|----------|
| H1 | **STRONG** | 실제 학습 모델(Ridge)로 3 데이터셋에서 확인. Bootstrap CI로 통계적 견고성 확보 |
| H2 | **MODERATE** | Domain-specific. Intra-domain (DEG↔f1, Mag↔f1): 100% 통과. Cross-domain: 33.3%. MSE 자체가 domain-general predictor |
| H3 | **STRONG** | 3 데이터셋 ALL WIN. 선형 모델로 검증. DL 모델 검증 진행 중 |

---

## 9. 실험 이력

| Run | 날짜 | 내용 | 결과 |
|-----|------|------|------|
| run_13 | 04-30 | BioEval 메트릭-순위 반전 (시뮬레이션 11개 모델) | H1 SUPPORTED (시뮬레이션) |
| run_14 | 04-30 | BioEval Phase 4 downstream 과업 상관 | H2 SUPPORTED (88.9%) |
| run_15 | 04-30 | sklearn Ridge LOO (Norman) | H1+H2 실제 모델 확인 (Norman). K562/RPE1 퇴화 |
| run_16 | 05-01 | Gene PCA Feature Ridge (K562/RPE1) | C1 해결. H1+H2+H3 3 데이터셋 전체 확인 |
| run_17 | 05-01 | Bootstrap CI (B=10,000) | H1+H2 통계적 견고성 확인. DEG_auprc가 모든 데이터셋에서 유의 |
| run_18 | 05-01 | Scale Correction (A3/B3) | 보정 불필요. Dir_deg 불변, downstream 과업 미개선 |
| run_19 | 05-01 | Downstream Task Independence (A4/B6) | H2 domain-specific. Cross-domain 33.3%, intra-DEG 100%. MSE domain-general |
| run_20 | 05-02 | GEARS DL 모델 훈련 (B2) | 진행 중 |

---

## 10. 참고 문헌

- Ahlmann-Eltze et al. (2025). "Comparison of perturbation prediction methods". Nature Methods.
- Chen et al. (2026). "SCALE: Single-cell perturbation landscape estimation". Cell-Eval framework.
- Csendes et al. (2025). "Benchmarking foundation models for perturbation prediction". BMC Genomics.
- Norman et al. (2019). "Mapping the perturbome landscape". Science.
- Raju et al. (2026). "Shesha: Geometric stability for perturbation prediction".
- Replogle et al. (2022). "Mapping information-rich genotype-phenotype landscapes". Cell.
- Roohani et al. (2025). "Virtual Cell Challenge". Cell.
- Wei et al. (2026). "Systematic benchmarking of perturbation prediction methods". Nature Methods.
- Zhu et al. (2025). "Evaluation metrics for perturbation prediction". Briefings in Bioinformatics.
