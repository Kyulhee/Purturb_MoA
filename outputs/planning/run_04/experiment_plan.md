# Run 04: Experiment Plan — 교세포 에피스태시스 전이 예측 + 소거실험

**Date:** 2026-04-29 | **Stage:** Planning (run_04) | **Framing:** outputs/framing/run_04/

---

## 1. 연구 질문 (from stages/02_framing.md run_04)

**"인과 불변성 기반 섭동 표현(z_tx)의 조합 예측 불확실성을 정량화할 때, (1) 불확실성이 조합 예측의 신뢰도를 유전자별·섭동별로 정량화할 수 있는가? (2) 한 세포유형의 에피스태시스 패턴과 불확실성으로 새 세포유형에서 에피스태시스가 발생할 조합을 예측할 수 있는가? (3) 불확실성 기반 능동학습이 교세포 에피스태시스 탐지 효율을 개선하는가?"**

### 하위 질문
1. **RQ1 (잔차 분해)**: 조합 섭동 예측 잔차를 에피스태시스 신호 vs 모델 오차 vs 노이즈로 원칙적 분해
2. **RQ2 (불확실성 정량화)**: MC Dropout 분산 + ICM 위반 점수를 결합하여 섭동 예측의 신뢰도를 유전자별·섭동별로 정량화
3. **RQ3 (교세포 에피스태시스 전이)**: 한 세포유형에서 관측한 에피스태시스 패턴으로, 새 세포유형에서 에피스태시스가 발생할 조합을 zero-shot으로 예측
4. **RQ4 (능동학습)**: 불확실성 기반 획득 함수가 교세포 에피스태시스 탐지에서 랜덤 선택 대비 효율 개선

### RQ3 변경 핵심 (vs run_03)
- run_03 RQ3: "잔차에서 에피스태시스 탐지" → AUROC=1.0 (동어반복: 잔차=에피스태시스)
- run_04 RQ3: "한 세포유형의 에피스태시스로 다른 세포유형의 에피스태시스 예측" → 순환 위험 제로

---

## 2. 방법론 아키텍처 (6-Phase)

```
Phase 1: FCR-ICM 모델 학습 (K562 + RPE1 각각)
  인코더 q(z|x,t) → z_x, z_t, z_tx 분해
  ICM 정규화: z_tx 세포 유형 불변성 강제 (MMD penalty)
  디코더 p(x|z_x, z_t, z_tx) → 재구성

Phase 2: 조합 예측 + 잔차 추출 (각 세포유형별)
  조합 함수: z_tx_A ⊕ z_tx_B → 예측 ŷ_AB (가법 기준)
  잔차: r = y_AB - ŷ_AB (유전자별 벡터)

Phase 3: 잔차 분해 (RQ1)
  r = r_epistasis + r_model + r_noise

Phase 4: 불확실성 정량화 (RQ2)
  U(g, comb) = w1·ICM_violation + w2·MC_dropout_var + w3·|r_additive|

Phase 5: 에피스태시스 점수 산출 (각 세포유형별, RQ3 입력)
  3-공식 민감도 분석 → 에피스태시스 점수 벡터 (K562용, RPE1용)
  ICM 가중 에피스태시스 점수

Phase 6: 교세포 에피스태시스 전이 예측 + 능동학습 (RQ3+RQ4)
  K562 에피스태시스 순위 → RPE1 에피스태시스 순위 예측
  ICM 정렬 z_tx 전이 vs 비인과 전이 vs 가법 잔차 전이
  불확실성 기반 능동학습 (교세포 설정)
```

---

## 3. 설계 결정과 근거

| 설계 결정 | 근거 | 대안 |
|-----------|------|------|
| FCR-ICM 기반 | z_x/z_t/z_tx 분해가 잔차 원천 구분에 구조적 이점 | CPA(잔차=오차), GEARS(GI 분류만) |
| 가법 조합 기준 | Pacalin 선례; trivial baseline과 명확한 분리 | 곱법(예측 기준은 가법이 표준) |
| 3-공식 민감도 분석 | Ajmal/Valenzuela 공식 격변 입증 | 단일 공식(재현성 위험) |
| ICM 위반 = 불확실성 신호 | 문헌에 없는 해석 (0개 경쟁자) | ICM 위반 = 오차(기존) |
| Replogle을 RQ3 주 데이터로 | K562+RPE1 공유 섭동 843개, 교세포 전이 평가 가능 | Norman(단일 세포유형, 전이 불가) |
| K562→RPE1 전이 방향 | K562 데이터가 더 큼(학습), RPE1이 독립 ground-truth | RPE1→K562(데이터 규모 불리) |
| 에피스태시스 순위 상관(rho)을 RQ3 주 지표로 | 연속 평가, 순환 불가, AUROC 동어반복 회피 | AUROC(순환 위험) |

---

## 4. 실험 설계

### 4.1 Phase 1: FCR-ICM 모델 학습

**아키텍처:**
- 인코더: [G] → MLP → z_x(16-d) + z_t(16-d) + z_tx(32-d)
- 디코더: [z_x, z_t, z_tx] → MLP → [G]
- G = 유전자 수 (Norman: ~2,000 HVG, Replogle: ~2,000 HVG)

**학습 (Replogle, RQ3 핵심):**
- K562와 RPE1 각각 독립 학습
- Loss = L_recon + β·KL + λ_ICM·MMD(z_tx^K562, z_tx^RPE1)
- ICM MMD: 두 세포유형의 z_tx 분포 정렬 (공유 843 섭동만 사용)

**학습 (Norman, RQ1 보조):**
- 단일 세포유형(K562)에서 학습
- ICM 정규화 없이 내부 일관성만 평가

**하이퍼파라미터:**

| 파라미터 | 범위 | 우선순위 |
|----------|------|----------|
| z_dim (z_tx) | [16, 32, 64] | 높음 |
| β (KL weight) | [0.1, 1.0] | 높음 |
| λ_ICM | [0.01, 0.1, 1.0, 10.0] | 높음 |
| learning rate | [1e-4, 3e-4] | 중간 |
| MC Dropout p | [0.05, 0.1, 0.2] | 높음 |

**데이터 분할:**
- Norman: 131 single KO (80/20 train/val), 104 double KO (test)
- Replogle: K562 843 섭동 (train), RPE1 843 섭동 (test, zero-shot)
  - 교세포 전이에서 RPE1은 완전 holdout (RQ3 순환 방지)

### 4.2 Phase 2: 조합 예측 + 잔차 추출

**Norman (단일 세포유형):**
- 가법 기준: ŷ_AB = ŷ_A + ŷ_B - μ
- 잔차: r_g = y_AB,g - ŷ_AB,g

**Replogle (교세포, RQ3 핵심):**
- 각 세포유형에서 독립적으로 조합 예측
- K562 조합 잔차 → K562 에피스태시스 점수
- RPE1 조합 잔차 → RPE1 에피스태시스 점수 (ground-truth)
- **주의**: Replogle은 단일 섭동만 있으므로, 조합 섭동은 가법 예측으로만 평가 가능
  - Norman의 104 double KO로 RQ1 검증
  - Replogle의 843 공유 단일 섭동으로 RQ3 에피스태시스 전이 평가

**RQ3 평가 설계 (핵심 변경):**

Replogle에는 조합 섭동 데이터가 없음. 에피스태시스 전이는 다음 방식으로 평가:

1. **유전자별 예측 오차 패턴 전이**: K562에서 유전자별 예측 오차가 큰 섭동 조합 순위 → RPE1에서도 동일 순서인지
2. **ICM 위반 패턴 전이**: K562에서 ICM 위반이 큰 섭동 → RPE1에서도 ICM 위반이 큰지
3. **잔차 구조 전이**: K562 잔차 공간 구조 → RPE1 잔차 공간 구조 상관

```
에피스태시스 전이 평가:
  source: K562의 843 섭동에 대한 에피스태시스 관련 점수
    - ICM 위반 점수 (z_tx 정렬 편차)
    - 예측 불확실성 (MC Dropout + ICM)
    - 가법 잔차 크기
  target: RPE1의 843 섭동에 대한 동일 점수
  평가: source 순위 vs target 순위의 Spearman rho
```

### 4.3 Phase 3: 잔차 분해 (RQ1)

**Norman 데이터에서 검증 (합성 + 실데이터):**

```
r_g = r_epistasis,g + r_model,g + r_noise,g
```

1. **r_noise,g**: VAE 재구성 분산 (MC Dropout 30회)
2. **r_model,g**: ICM 위반 기여 (ICM_violation × ∂ŷ_g/∂z_tx)
3. **r_epistasis,g**: r_g - r_model,g - r_noise,g

**합성 검증 (run_09에서 r=1.000 달성, 재확인):**
- 알려진 에피스태시스 구조 주입 → 분해 정확도 측정

**Norman 실데이터 간접 검증:**
- GEARS GI 분류와 분해 결과의 일치도
- 알려진 상호작용(DUSP9+ETS2 등) 탐지 여부

### 4.4 Phase 4: 불확실성 정량화 (RQ2)

**다중 신호원 결합:**

```
U(g, comb) = w1 · σ_ICM(g, comb) + w2 · σ_MC(g, comb) + w3 · |r_add(g, comb)|
```

**Replogle 교세포 UQ (RQ2+RQ3 연결):**
- K562에서 UQ 산출 → RPE1에서 UQ 산출
- UQ 순위 전이: ρ(U^K562, U^RPE1) — 불확실성도 세포유형 간 보존되는가?

**평가:**

| 지표 | Baseline | 타겟 | 데이터 |
|------|----------|------|--------|
| U-Error Spearman rho | ~0.3-0.5 | > 0.6 | Replogle |
| Coverage (90% CI) | — | 0.85-0.95 | Replogle |
| UQ 전이 rho | — | > 0.4 (보조) | Replogle K562→RPE1 |

### 4.5 Phase 5: 에피스태시스 점수 산출 (RQ3 입력)

**3-공식 민감도 분석 (Norman + Replogle):**

| 공식 | Expected (null) | 에피스태시스 정의 |
|------|----------------|-------------------|
| Additive | ŷ_A + ŷ_B - μ | r = observed - expected |
| Multiplicative | ŷ_A × ŷ_B / μ | r = observed/expected - 1 |
| Product neutrality | fitness_A × fitness_B | r = fitness_AB/(fitness_A·fitness_B) - 1 |

**에피스태시스 점수 벡터 (각 섭동별):**
- E_add(s), E_mult(s), E_prod(s) — 3공식 각각의 에피스태시스 점수
- ICM 가중: E_ICM(s) = ICM_weight(s) × mean(E_add, E_mult, E_prod)
- Formula agreement: 3공식 일치율 (run_09: 76.6%)

**교세포 전이 (RQ3 핵심):**

```
K562 에피스태시스 점수 순위 → RPE1 에피스태시스 점수 순위 예측

전이 방법:
1. ICM 정렬 전이 (본 제안): ICM 정규화된 z_tx 기반 에피스태시스 순위 전이
2. 비인과 전이: ICM 없이 z_tx 기반 전이 (음의 상관 예상, run_07: -0.35)
3. 가법 잔차 순위 전이: 단순 잔차 크기 순위 전이
4. UQ 순위 전이: 불확실성 순위 전이
5. 랜덤: 우연 수준

핵심 비교: ICM 전이 vs 비인과 전이 → ICM 정규화의 전이 개선율
```

### 4.6 Phase 6: 능동학습 (RQ4)

**교세포 능동학습 (새로운 설정):**

```
시나리오: K562에서 모든 섭동을 관측, RPE1에서 k개만 관측 가능
목표: 불확실성 기반으로 RPE1에서 관측할 섭동 선택 → 에피스태시스 탐지 효율 극대화

획득 함수:
  α(s) = U_K562(s) + λ · transfer_gain(s)

  - U_K562(s): K562에서의 불확실성 (전이 가능한 불확실성)
  - transfer_gain(s): ICM 정렬도 × 에피스태시스 예상치

비교:
  1. 불확실성 + ICM 전이 (본 제안)
  2. 불확실성만 (no ICM)
  3. 랜덤 선택
  4. 엔트로피 기반 선택
```

**평가:**

| 지표 | Baseline | 타겟 | 데이터 |
|------|----------|------|--------|
| AL Top-k 개선율 (교세포) | k/n (random) | > 2× random | Replogle |
| AL Transfer overlap | — | > 0.5 | Replogle |
| AL Top-k recall (Norman) | k/n | > 2× random | Norman |

---

## 5. 소거 실험 매트릭스 (핵심: ICM 전이 소거)

| ID | FCR | ICM | 3-공식 | 잔차 분해 | UQ 결합 | 전이 방법 | 예상 결과 |
|----|-----|-----|--------|----------|---------|----------|----------|
| A1 | ✓ | ✓ | ✓ | ✓ | ✓ | ICM 전이 | Full model |
| A2 | ✓ | — | ✓ | ✓ | ✓ | 비인과 전이 | ICM 정규화 기여 측정 (**핵심**) |
| A3 | ✓ | ✓ | — (단일) | ✓ | ✓ | ICM 전이 | 공식 민감도 기여 |
| A4 | ✓ | ✓ | ✓ | — (trivial) | ✓ | ICM 전이 | 잔차 분해 기여 |
| A5 | ✓ | ✓ | ✓ | ✓ | — (MC only) | ICM 전이 | 다중 신호원 기여 |
| A6 | ✓ | ✓ | ✓ | ✓ | ✓ | — (no AL) | AL 기여 |
| A7 | ✓ | ✓ | ✓ | ✓ | ✓ | 잔차 순위 전이 | 전이 방법 기여 (**핵심**) |
| A8 | ✓ | ✓ | ✓ | ✓ | ✓ | UQ 순위 전이 | UQ 전이 기여 |
| B1 | CPA | — | — | — | MC only | — | CPA baseline |
| B2 | GEARS | — | — | — | — | — | GEARS baseline |

**핵심 소거 (A1 vs A2, A1 vs A7):**
- A1 vs A2: ICM 정규화가 에피스태시스 전이를 개선하는가? → **ICM 전이 개선율**
- A1 vs A7: ICM 전이가 가법 잔차 순위 전이보다 우수한가? → **방법론 기여**

---

## 6. 평가 지표 요약

| RQ | 지표 | Baseline | 타겟 | 데이터 | 순환 위험 |
|----|------|----------|------|--------|----------|
| RQ1 | 잔차 분해 정확도 (합성) | — | r > 0.7 | 합성 | 없음 |
| RQ2 | U-Error Spearman rho | ~0.3-0.5 | > 0.6 | Replogle | 없음 |
| RQ2 | Coverage (90% CI) | — | 0.85-0.95 | Replogle | 없음 |
| RQ3 | Cross-CT 에피스태시스 rho | 0 (우연) | > 0.4 | Replogle K562→RPE1 | **없음** |
| RQ3 | Top-k overlap (에피스태시스) | k/n (우연) | > 2× random | Replogle | 없음 |
| RQ3 | ICM 전이 개선율 | Baseline(no ICM) rho | > 1.5× baseline | Replogle | 없음 |
| RQ3 | 에피스태시스 Precision | ~0.3 (trivial) | > 0.6 | Norman | 낮음 |
| RQ4 | AL Top-k 개선율 (교세포) | k/n (random) | > 2× random | Replogle | 없음 |
| RQ4 | AL Transfer overlap | — | > 0.5 | Replogle | 없음 |

---

## 7. 베이스라인 계층 (RQ3 특화)

| 계층 | 방법 | 구현 | 예상 성능 | 순환 위험 |
|------|------|------|----------|----------|
| 우연 | 랜덤 조합 순서 | rho ≈ 0, overlap ≈ k/n | — | 없음 |
| 비인과 전이 | ICM 없이 z_tx로 전이 | A2 소거 | rho < 0 (run_07: -0.35) | 없음 |
| 가법 잔차 전이 | 단순 잔차 크기 순위 전이 | A7 소거 | rho > 0 (부분 상관) | 없음 |
| UQ 순위 전이 | 불확실성 순위 전이 | A8 소거 | rho > 0 | 없음 |
| **ICM 전이 (본 제안)** | ICM 정렬 z_tx 기반 | A1 | rho > 0.4, > 1.5× baseline | 없음 |

---

## 8. 데이터

| 데이터셋 | 용도 | 접근성 | 전처리 |
|----------|------|--------|--------|
| Norman (2019) | RQ1 잔차 분해 검증 | GEO GSE133344 | Scanpy, 2000 HVG |
| Replogle (2022) | RQ2-4 주 검증 | GEO GSE142398 | K562 + RPE1 분리, 843 공유 섭동 |
| PORTAL (2026) | RQ3 대규모 외부 검증 (선택) | 접근성 확인 필요 | — |

**Replogle 전처리 (RQ3 핵심):**
1. K562와 RPE1 각각 독립 전처리
2. 공유 843 섭동 식별
3. 각 세포유형에서 2000 HVG 선택 (공통 유전자 사용)
4. Control 셀 분리
5. 정규화 + 로그 변환

---

## 9. 핵심 리스크와 완화

| 리스크 | 확률 | 영향 | 완화 |
|--------|------|------|------|
| 에피스태시스 교세포 보존이 우연 수준 (rho < 0.2) | 중간 | 높음 | rho 0.2-0.4는 partial success로 처리; ICM 전이 개선율이 1.5× 이상이면 moderate claim 가능 |
| ICM 전이가 비인과 전이보다 나쁨 | 낮음 | 치명적 | A1 vs A2 소거로 명확히 측정; 음의 상관이면 ICM 정렬이 해로움을 보고 |
| Precision이 trivial baseline과 유의미 차이 없음 | 중간 | 높음 | Precision만으로 판단하지 않고 rho + overlap + 개선율 종합 평가 |
| Replogle에 조합 섭동 없음 | 확실 | 중간 | 단일 섭동 기반 에피스태시스 프록시(예측 오차, ICM 위반)로 전이 평가; Norman에서 조합 섭동 검증 |
| 공식 민감도가 교세포 전이에 치명적 | 낮음 | 중간 | 3-공식 각각 교세포 전이 rho 측정; 불일치 자체를 결과로 보고 |
| 단일 아키텍처 한계 (FCR-ICM만) | 높음 | 중간 | B1(CPA), B2(GEARS) baseline 비교; BuDDI/scDRP 비교는 선택 |

---

## 10. 이전 run 실패에서의 설계 원칙

| 교훈 | 출처 | 본 설계 반영 |
|------|------|-------------|
| R2만으로 novelty 불충분 | run_02/08 | Precision/Recall/rho/overlap로 평가 전환 |
| 강한 trivial baseline 주의 | run_08/09 | Precision 측정 포함; trivial baseline 명시 |
| AUROC=1.0 동어반복 | run_09 | RQ3를 교세포 전이로 재정의 (순환 불가) |
| 소거실험 미실행 | run_09 | A2-A8 전체 소거 매트릭스 설계 |
| Precision 미측정 | run_09 | RQ3에 Precision 지표 추가 |
| ICM 전이 R2는 높아도 novelty 약함 | run_08 | R2 대신 rho + 개선율로 평가 |
| Comp loss 불필요 | run_06 | 가법 기준, 복잡한 loss 불필요 |
| 잠재공간-유전자공간 갭 | run_05 | 모든 평가 유전자 공간에서 수행 |

---

## 11. 컴퓨팅 자원 추정

| 단계 | 예상 시간 | 비고 |
|------|----------|------|
| Norman + Replogle 전처리 | 2-3h | Scanpy, 메모리 16GB+ |
| FCR-ICM 학습 (2 세포유형) | 4-8h | GPU 필수, RTX 4060 Ti |
| 잔차 분해 + UQ | 2-3h | MC Dropout 30회 |
| 3-공식 에피스태시스 산출 | 1h | 공식별 계산 |
| 교세포 전이 평가 | 1-2h | 순위 상관, overlap |
| AL 시뮬레이션 | 3-5h | 교세포 + Norman |
| 소거실험 (10 구성) | 2-3일 | × 3-fold CV |

**총 예상: 4-6일 (GPU 포함)**

---

## 12. 논문 구조 제안

1. **Introduction**: Perturb-seq 조합 예측 실패 → 잔차에서 무엇을 배울 수 있는가? 에피스태시스가 세포유형 간에 보존되는가?
2. **Related Work**: 섭동 예측, 에피스태시스 정의, 인과 표현 학습(FCR/ICM/IEM), 불확실성(CIPHER)
3. **Methods**: FCR-ICM 잔차 분해 → 다중 신호원 UQ → 교세포 에피스태시스 전이 → 능동학습
4. **Results**:
   - 4.1 RQ1: 잔차 분해가 에피스태시스 신호를 정밀하게 분리하는가?
   - 4.2 RQ2: 다중 신호원 불확실성이 예측 오차와 상관하는가?
   - 4.3 RQ3: 한 세포유형의 에피스태시스로 다른 세포유형의 에피스태시스를 예측할 수 있는가? (ICM 전이 개선율)
   - 4.4 RQ4: 불확실성 기반 교세포 능동학습의 효율
5. **Discussion**: 에피스태시스 교세포 보존의 의미, ICM 불변성의 확장, 한계(단일 아키텍처, moderate rho), 확장(PORTAL)
