# 잔차→에피스태시스 분해: 이론적 프레임워크 재정립

**Date:** 2026-04-29 | **Context:** novelty assessment에서 식별된 선행 논문들에 대한 심층 비교

---

## 1. 식별된 선행 논문 (분야 간 교차 검증)

### 1.1 직접적 선행 — 생물학

| 논문 | 분야 | 핵심 방법 | 우리와의 관계 |
|------|------|----------|-------------|
| **Pacalin et al.** (2025, Nature Biotech) | CRISPRai 이중섭동 | 단일 섭동 log2FC 가법 합 → 잔차로 GI 분류 | **최근접** — 단 baseline이 훈련된 모델 아님 |
| **MoCHI** (Faure & Lehner 2024, Genome Biology) | 단백질 DMS | NN으로 interpretable model 적합, 잔차에서 energetic coupling/epistasis 추정 | **유사** — 모델 적합 후 잔차 분석. 단 DMS 도메인, biophysical model 기반 |
| **Otwinowski et al.** (2018, PNAS) | 단백질 DMS | I-spline 비선형 매핑 → 잠재 형질 가법 모델 → HOC 잔차 = 에피스태시스 | **유사** — 가법 모델 잔차→고차 결합. 단 단백질, 잠재형질 공간 |
| **Birgy et al.** (2026, Nature Comms) | 단백질 안정성 | 안정성 모델 편차 = 구조적 에피스태시스 | **부분 유사** — 모델 잔차를 에피스태시스로 해석. 단 3D 구조 기반 |

### 1.2 간접적 선행 — ML/XAI

| 논문 | 분야 | 핵심 방법 | 우리와의 관계 |
|------|------|----------|-------------|
| **Diamond** (Chen et al. 2024) | ML 해석가능성 | ML 모델에서 비가법적 feature interaction 발견, knockoff로 FDR 제어 | **개념적 유사** — ML 모델의 비가법적 효과 탐지. 단 feature interaction, 생물학적 에피스태시스 아님 |
| **Tan et al.** (2018, JMLR) | XAI 이론 | 가법 설명이 비가법 모델에서 무엇을 놓치는지 분석. 잔차=비가법적 구조 | **이론적 선행** — "가법 근사의 잔차 = 비가법적 신호"라는 우리 핵심 아이디어의 정확한 이론적 기반 |
| **"Do Not Trust Additive Explanations"** (Gosiewska & Biecek 2019) | XAI | 가법 설명의 불신, 비가법적 효과의 중요성 입증 | **개념적 지지** — 가법 근사의 한계를 보이는 것은 우리의 동기와 같음 |
| **RED** (Brugger et al. 2025) | 방정식 발견 | 모델 잔차를 기반으로 방정식 개선. 잔차에 구조가 있으면 새 항 발견 | **방법론적 유사** — 잔차에 구조가 있으면 그것은 "신호"임. 단 방정식 발견 도메인 |

### 1.3 간접적 선행 — 인과추론/전이

| 논문 | 분야 | 핵심 방법 | 우리와의 관계 |
|------|------|----------|-------------|
| **C3TL** (arXiv:2603.13051) | 인과 전이 | 인과 불변성으로 cross-domain 전이 | **ICM 선행** — 같은 원리, 다른 도메인 (bulk) |
| **Subbaswamy et al.** (2019) | dataset shift | 불변 예측 분포로 외부 타당성 보장 | **이론적 기반** — ICM의 건전성 뒷받침 |

---

## 2. 핵심 질문: 우리의 novelty는 어디에 있는가?

### 2.1 아이디어 "가법 모델 잔차 = 비가법적 신호"의 선행성 검토

**이 아이디어는 새로운가?** — **아니다.** 

회귀분석에서 교호작용항(interaction term)을 잔차에서 발견하는 것은 **통계학의 기본 교리**다. ANOVA에서 main effect 잔차에서 interaction effect를 검출하는 것이 표준이다.

Tan et al. (2018)이 이를 XAI 맥락에서 공식화: "가법 설명은 비가법적 모델의 잔차를 놓친다." Diamond (2024)가 이를 ML 모델에 적용: "ML 모델에서 비가법적 feature interaction을 FDR 제어와 함께 발견."

**우리의 차별화는 아이디어 자체가 아니라:**

1. **적용 도메인**: 섭동 생물학에서 조합 예측 잔차 → 에피스태시스 (이 도메인에서 0개)
2. **분해 구조**: 단순 "잔차=에피스태시스"가 아닌, FCR의 z_x/z_t/z_tx 분해가 잔차의 **원천 식별**을 가능하게 함
3. **3-공식 민감도**: 에피스태시스 정의 자체가 공식에 의존한다는 것을 명시적으로 다룸
4. **ICM과의 시너지**: ICM 정규화가 잔차 분해를 가능하게 하는 전제조건 (z_tx 불변성 없이는 교세포 혼란이 잔차에 묻힘)

### 2.2 재정립된 novelty 계층도

```
Layer 3 (우리만): 섭동 예측 잔차 → 에피스태시스 분해 + ICM + 3-공식 + AL
                  ↑ 적용 도메인 + 통합 프레임워크
Layer 2 (선행 있음): 가법 모델 잔차 → 비가법적 신호 탐지
                      ↑ Tan 2018, Diamond 2024, Otwinowski 2018, MoCHI 2024
Layer 1 (일반 원리): 모델 오류 = 규격 오차 + 누락된 구조 + 노이즈
                      ↑ 회귀분석, ANOVA, 통계학의 기본
```

---

## 3. 재정립된 이론적 프레임워크

### 3.1 공식화

FCR-ICM의 조합 예측에서:

```
y_obs = y_pred + r
r = r_epistasis + r_model + r_noise
```

여기서:
- **r_epistasis**: 가법적 조합 가정에 의해 포착되지 않은 **진짜 비가법적 생물학적 효과**
- **r_model**: 모델의 규격 오차 (ICM 불충분, 디코더 용량 부족, 학습 불충분)
- **r_noise**: 측정/생물학적 노이즈

이 분해의 **핵심 착상**: FCR의 구조적 분해(z_x/z_t/z_tx)가 없으면 r_epistasis와 r_model을 구별할 방법이 없다.

**이론적 근거:**

1. **Tan et al. (2018)의 정리**: 가법 근사 g(x)가 비가법적 함수 f(x)를 근사할 때, 잔차 f(x)-g(x)는 비가법적 구조를 포함한다. 단, 이 잔차는 "누락된 구조"이지 반드시 "에피스태시스"는 아니다.

2. **우리의 확장**: FCR 구조가 잔차의 원천을 구별:
   - r_model의 하한은 ICM violation score로 추정 (z_tx가 cell type에 의존적이면 모델이 불충분)
   - r_model의 상한은 OOD distance로 추정 (학습 분포에서 먼 조합일수록 모델 오차 가능성 증가)
   - r_epistasis는 |r| - r_model_bound로 추정

3. **Valenzuela (2025)의 공식 민감도**: r_epistasis의 부호와 크기가 공식(additive/multiplicative/product)에 의존. → 3-공식 일치도(0.766)가 에피스태시스 신호의 **견고성** 지표.

### 3.2 Diamond와의 비교 — 핵심 차이

| 측면 | Diamond (2024) | FCR-ICM (우리) |
|------|---------------|---------------|
| 목적 | ML 모델에서 feature interaction 발견 | 섭동 예측에서 에피스태시스 분해 |
| 대상 | 일반 ML (tabular, DNN, transformer) | 생물학적 섭동 예측 모델 |
| 모델 | 모델 불가지 (black-box) | 구조적 모델 (FCR 분해) |
| 잔차 분해 | ❌ — interaction 있/없음만 판정 | ✅ — epistasis/model error/noise 3분해 |
| FDR 제어 | ✅ — knockoff 기반 | ❌ — 현재 없음 (약점!) |
| 비가법성 정의 | feature interaction | 3-공식 에피스태시스 (add/mult/prod) |
| 불확실성 | FDR만 | OOD + MC Dropout + ICM |
| 능동학습 | ❌ | ✅ |

**핵심 차이**: Diamond는 "interaction이 있는가?"를 FDR 제어와 함께 대답. 우리는 "잔차를 무엇으로 구성되어 있는가?"를 분해. **질문이 다름.**

### 3.3 MoCHI와의 비교 — 핵심 차이

| 측면 | MoCHI (2024) | FCR-ICM (우리) |
|------|-------------|---------------|
| 도메인 | 단백질 DMS | 유전자 섭동 (Perturb-seq) |
| 모델 | interpretable biophysical model | VAE + ICM |
| 잔차 분석 | 잔차에서 global nonlinearity 추정 | 잔차를 3성분으로 분해 |
| 고차 에피스태시스 | ✅ (higher-order terms 지원) | ❌ (2차 조합만) |
| 교세포 전이 | ❌ | ✅ (ICM) |
| 표현 공간 | 에너지 (물리적) | z_tx (학습된) |
| 검증 | biophysical ground truth 있음 | 실데이터 ground truth 없음 (약점!) |

**핵심 차이**: MoCHI는 단백질 물리 모델 기반으로 ground truth가 있음. 우리는 데이터 기반 학습 표현이므로 직접 검증 불가. 대신 ICM+OOD로 간접 검증.

---

## 4. Novelty 재평가

### 수정된 novelty (이론적 선행 반영)

1. **잔차→에피스태시스 분해** — ★★★★★ → ★★★★☆
   - 아이디어(가법 잔차=비가법 신호)는 통계학의 기본이고 Tan/Diamond가 ML에 적용
   - **우리만의 것**: 섭동 생물학 적용 + FCR 구조적 분해로 원천 식별 + 3-공식 민감도
   - **약화 요인**: Diamond가 "ML 모델에서 비가법적 효과 탐지"를 이미 수행

2. **ICM 교세포 전이** — ★★★☆☆ (변화 없음, C3TL 선행)

3. **3-공식 에피스태시스 민감도** — ★★★★☆ → ★★★☆☆
   - Valenzuela가 product neutrality 이론 제시
   - Ajmal/Chitra가 공식 민감도 입증
   - **우리만의 것**: 예측 잔차에 3-공식을 적용하고 일치도를 정량화
   - **약화 요인**: 공식 민감도 자체는 이미 알려진 사실

4. **OOD distance 기반 UQ** — ★★★☆☆ (변화 없음)

5. **능동학습** — ★★★☆☆ (변화 없음)

### 수정된 통합 novelty

| 기준 | 평가 |
|------|------|
| 아이디어 독창성 | ★★★☆☆ — "잔차→비가법적 신호"는 일반적 원리 |
| 적용 독창성 | ★★★★★ — 섭동 생물학에서 0개 |
| 방법론 독창성 | ★★★★☆ — FCR 분해로 원천 식별은 새로움 |
| 통합 독창성 | ★★★★☆ — 6개 기능을 하나의 프레임워크에 통합 |

---

## 5. 이론적 위치 재설정

### 논문에서 써야 할 포지셔닝

**올바른 포지셔닝:**
> "가법 모델의 잔차에 구조가 있다는 것은 새로운 발견이 아니다 (Tan 2018, Diamond 2024). 우리의 기여는 **섭동 생물학**에서 이 원리를 적용할 때, (1) 잔차의 원천을 FCR 분해로 식별하고, (2) 에피스태시스 정의의 공식 의존성을 명시적으로 다루며, (3) ICM 정규화가 분해 가능성의 전제조건임을 보이는 것이다."

**피해야 할 포지셔닝:**
> "예측 잔차에서 에피스태시스를 분해하는 것은 우리가 최초다" — 정확하지 않음. Otwinowski, MoCHI, Pacalin이 유사 접근. Diamond가 ML에서 일반화.

### 논문 Introduction에 들어가야 할 선행 인용

1. **Tan et al. (2018)**: "가법 근사의 잔차는 비가법적 구조를 포함한다" — 우리 이론의 일반적 기반
2. **Diamond (Chen et al. 2024)**: "ML 모델에서 비가법적 feature interaction을 FDR 제어와 함께 발견" — ML 선행
3. **MoCHI (Faure & Lehner 2024)**: "DMS에서 잔차에서 energetic coupling 추정" — 생물학 선행
4. **Pacalin et al. (2025)**: "CRISPRai에서 가법 합 잔차로 GI 분류" — 최근접 섭동 생물학 선행
5. **Valenzuela (2025)**: "Product neutrality의 기계적 근거" — 에피스태시스 공식 이론

### 우리의 차별화 진술 (Differentiation Statement)

> Diamond는 black-box ML 모델에서 feature interaction의 존재 여부를 FDR 제어와 함께 판정하지만, 잔차를 구조적으로 분해하지는 않는다. MoCHI는 DMS 잔차에서 energetic coupling을 추정하지만, biophysical model 기반이므로 model error와 epistasis의 구분이 불필요하다. Pacalin은 가법 합 잔차로 GI를 분류하지만, 훈련된 예측기가 아닌 단순 합의 잔차이므로 model error가 없다.
>
> **우리의 기여**: 훈련된 신경망 예측기의 잔차에서 **에피스태시스와 모델 오차를 구별**하는 원칙적 방법을 제안한다. 이 구별은 FCR의 표현 분해(z_x/z_t/z_tx)와 ICM 정규화가 가능하게 하는 것으로, 구조적 모델 없이는 불가능하다.

---

## 6. 약점 보완 방안

### 6.1 FDR 제어 (Diamond 대비 약점)
- Diamond의 knockoff 프레임워크를 우리 잔차 분해에 적용 가능한지 탐색
- 최소한: 순열 검정으로 에피스태시스 검출의 FDR 추정

### 6.2 실데이터 검증 (MoCHI 대비 약점)
- PORTAL 665K 쌍을 holdout 검증으로 사용
- 단기적으로: Norman 데이터에서 cross-validation 기반 간접 검증

### 6.3 고차 에피스태시스 (MoCHI 대비 약점)
- 현재 2차 조합만 지원
- 장기적으로: triple-KO 등 고차 조합으로 확장

---

## 7. 결론: 수정된 논문 가능성

**옵션 A(잔차 분해 + ICM 전이) 수정 가능성: 30-40%** (기존 35-45%에서 하향)

하향 이유:
- "잔차→비가법 신호" 아이디어가 Tan(2018)/Diamond(2024)에서 이미 일반화됨
- MoCHI/Otwinowski가 생물학(단백질)에서 이미 유사 접근
- 우리의 차별화가 "적용 도메인"과 "FCR 구조적 분해"에 있으나, 이것만으로는 top venue에 부족할 수 있음

**논문 가능하게 만드는 핵심 보완:**
1. **FDR 제어 추가** — Diamond의 knockoff를 잔차 분해에 적용. 이렇게 되면 "에피스태시스 발견에 FDR 보장"이라는 강력한 차별화
2. **Replogle 다세포유형에서 일관성 검증** — ICM 진가 + RQ2 개선 동시 달성
3. **PORTAL 대규모 검증** — 665K 쌍으로 precision/recall 정량화
