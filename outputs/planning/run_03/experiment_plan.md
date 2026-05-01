# Run 03: Experiment Plan — 에피스타시스 탐지 + 불확실성 정량화

**Date:** 2026-04-28 | **Stage:** Planning (run_03) | **Framing:** outputs/framing/run_03/

---

## 1. 연구 질문 (from stages/02_framing.md)

**"인과 분해된 섭동 표현(z_tx)의 조합 예측 불확실성을 정량화할 때, (1) 예측 잔차에서 비가법적 유전자 상호작용(에피스타시스)을 체계적으로 분리할 수 있는가? (2) 이 불확실성이 제한된 실험 예산 하에서 어떤 조합 섭동을 우선 실험해야 할지 안내할 수 있는가?"**

### 하위 질문
1. **RQ1 (잔차 분해)**: 조합 섭동 예측 잔차를 에피스타시스 신호 vs 모델 오차 vs 노이즈로 원칙적으로 분해
2. **RQ2 (불확실성 정량화)**: ICM 위반 점수, MC Dropout 분산, 가법 잔차를 결합한 유전자별·섭동별 불확실성
3. **RQ3 (에피스타시스 탐지)**: 분해된 잔차에서 에피스타시스 존재·강도·유형 탐지
4. **RQ4 (능동학습)**: 불확실성 기반 획득 함수의 조합 섭동 실험 효율 개선

---

## 2. 방법론 아키텍처

```
Phase 1: FCR-ICM 기반 모델 학습
  인코더 q(z|x,t) → z_x, z_t, z_tx 분해
  ICM 정규화: z_tx 세포 유형 불변성 강제 (MMD penalty)
  디코더 p(x|z_x, z_t, z_tx) → 재구성

Phase 2: 조합 예측 + 잔차 추출
  조합 함수: z_tx_A ⊕ z_tx_B → 예측 ŷ_AB (가법 기준)
  잔차: r = y_AB - ŷ_AB (유전자별 벡터)

Phase 3: 잔차 분해 (RQ1 핵심 기여)
  r = r_epistasis + r_model + r_noise
  - r_epistasis: 비가법적 상호작용 신호 (3공식 민감도 분석)
  - r_model: ICM 위반 점수 × 가중치 (모델 오지정)
  - r_noise: VAE 재구성 분산 (생물학적/기술적 노이즈)

Phase 4: 불확실성 정량화 (RQ2)
  U(g, comb) = w1·ICM_violation + w2·MC_dropout_var + w3·|r_additive|
  - 유전자별(g), 섭동별(comb) 분해 가능
  - 보정: conformal calibration on held-out combinations

Phase 5: 에피스타시스 탐지 (RQ3)
  3-공식 민감도 분석:
    expected_add = ŷ_A + ŷ_B - 1  (가법)
    expected_mult = ŷ_A × ŷ_B     (곱법)
    expected_prod = fitness_A × fitness_B  (Product neutrality)
  공식 간 불일치 → "epistemic uncertainty from formula choice"
  3공식 모두 일치 → "high-confidence epistasis"

Phase 6: 능동학습 시뮬레이션 (RQ4)
  획득 함수: α(comb) = U(comb) + λ·diversity(comb)
  Norman 데이터에서 순차적 조합 선택 시뮬레이션
```

---

## 3. 설계 결정과 근거

| 설계 결정 | 근거 | 대안 검토 |
|-----------|------|----------|
| FCR-ICM을 기반 모델로 채택 | z_x/z_t/z_tx 분해가 잔차 원천 구분에 구조적 이점 | CPA(잔차=오차 처리), GEARS(GI 분류만) |
| 가법 조합을 기준 예측으로 사용 | Pacalin 2025가 가법 합으로 GI 분류한 선례; trivial baseline과의 명확한 분리 | 곱법 조합(Valenzuela 근거이나 예측 기준은 가법이 표준) |
| 3-공식 민감도 분석 포함 | Chitra 2025, Ajmal 2025가 공식 선택에 따른 결과 격변 입증; Valenzuela 2025가 Product neutrality 이론적 우수성 입증 | 단일 공식 사용(재현성 위험) |
| ICM 위반을 불확실성 신호로 해석 | ICM 위반 = z_tx가 문맥 의존적 = 잠재적 에피스타시스 신호; 이 해석은 문헌에 없음 (0개 경쟁자) | ICM 위반 = 버려야 할 오차(기존 해석) |
| Norman을 주 검증 데이터로 사용 | 131 single + 104 double KO, ground-truth 조합 데이터; GEARS 등 기존 방법의 표준 벤치마크 | Replogle(조합 섭동 없음), PORTAL(규모는 크나 접근성 미확) |
| MC Dropout + ICM 위반 + 가법 잔차 다중 신호원 | 단일 불확실성 신호로는 한계; CIPHER도 Bayesian 단일 신호원 | Deep Ensemble(계산 비용 5×), Conformal alone(적응성 부족) |

---

## 4. 실험 설계

### 4.1 Phase 1: FCR-ICM 모델 학습

**아키텍처:**
- 인코더: [G] → MLP → z_x(16-d) + z_t(16-d) + z_tx(32-d)
- 디코더: [z_x, z_t, z_tx] → MLP → [G]
- G = 유전자 수 (Norman: ~2,000 HVG)

**학습:**
- Loss = L_recon + β·KL + λ_ICM·MMD(z_tx^A, z_tx^B)
- L_recon: MSE (유전자 공간)
- KL: β-VAE 항 (β ∈ {0.1, 1.0})
- MMD: RBF 커널, multi-scale (γ ∈ {0.1, 1.0, 10.0})

**하이퍼파라미터 탐색:**

| 파라미터 | 범위 | 우선순위 |
|----------|------|----------|
| z_dim (z_tx) | [16, 32, 64] | 높음 |
| β (KL weight) | [0.1, 1.0] | 높음 |
| λ_ICM | [0.01, 0.1, 1.0, 10.0] | 높음 |
| learning rate | [1e-4, 3e-4] | 중간 |
| encoder layers | [2, 3] | 낮음 |
| hidden dim | [128, 256] | 낮음 |

**데이터 분할:**
- Norman 데이터: 131 single KO + 104 double KO
- Single KO: 80% train / 20% val
- Double KO: **전체를 test에 사용** (조합 예측은 기본적으로 OOD)
- 교차 검증: 5-fold on single KO (double KO는 고정 test)

### 4.2 Phase 2: 조합 예측 + 잔차 추출

**조합 함수:**
- 기준: 가법 ŷ_AB = μ + (ŷ_A - μ) + (ŷ_B - μ) = ŷ_A + ŷ_B - μ
  - μ = control 평균, ŷ_A = 단일 섭동 A 예측, ŷ_B = 단일 섭동 B 예측
- z_tx 공간: z_tx_AB = z_tx_A + z_tx_B (가법 조합)
- 디코더로 유전자 공간 예측 복원

**잔차:**
- r_g = y_AB,g - ŷ_AB,g (유전자 g에 대한 잔차)
- 정규화: r_g / σ_g (유전자별 분산으로 스케일)

### 4.3 Phase 3: 잔차 분해 (RQ1)

**분해 방법:**

```
r_g = r_epistasis,g + r_model,g + r_noise,g
```

1. **r_noise,g**: VAE 재구성 분산
   - MC Dropout 30회 → Var(ŷ) per gene
   - 분산이 큰 유전자 = 노이즈 높은 유전자

2. **r_model,g**: ICM 위반 기여
   - ICM_violation_AB = MMD(z_tx_AB^train, z_tx_AB^test) 또는
   - ICM_violation_AB = ||E[z_tx|A, train] - E[z_tx|A, test]||²
   - r_model,g ∝ ICM_violation_AB × ∂ŷ_g/∂z_tx

3. **r_epistasis,g**: 잔차 - 모델 오차 - 노이즈
   - r_epistasis,g = r_g - r_model,g - r_noise,g
   - 이 성분이 비가법적 상호작용 신호

**합성 검증 (RQ1 타겟: 분해 성분과 ground-truth r > 0.7):**
- 합성 데이터: 알려진 에피스타시스 구조 주입
  - y = f(x) + ε_epistasis + ε_model + ε_noise (각 성분 독립)
  - 분해 정확도: Corr(r̂_epistasis, ε_epistasis) 등
- Norman 실데이터: 간접 검증
  - GEARS GI 분류와 분해 결과의 일치도
  - 알려진 상호작용(DUSP9+ETS2 시너지 등) 탐지 여부

### 4.4 Phase 4: 불확실성 정량화 (RQ2)

**다중 신호원 결합:**

```
U(g, comb) = w1 · σ_ICM(g, comb) + w2 · σ_MC(g, comb) + w3 · |r_add(g, comb)|
```

- σ_ICM: ICM 위반 점수의 유전자별 투영
- σ_MC: MC Dropout 예측 분산 (30회 forward pass)
- |r_add|: 가법 잔차의 절댓값
- w1, w2, w3: 검증 세트에서 최적화 (grid search)

**보정:**
- Conformal prediction on held-out combinations
- 90% CI coverage 타겟: 0.85-0.95
- 분위수 회귀 또는 split-conformal 적용

**평가:**
- Uncertainty-Error Spearman rho: ρ(U, |error|) — 불확실성이 높은 곳이 실제 오차도 큰가?
- Coverage: 90% CI가 실제로 85-95% 포함하는가?
- 유전자별 구분력: 예측 가능한 유전자 vs 불가능한 유전자의 U 분리도

### 4.5 Phase 5: 에피스타시스 탐지 (RQ3)

**3-공식 민감도 분석:**

| 공식 | Expected (null) | 에피스타시스 정의 | 근거 |
|------|----------------|-------------------|------|
| Additive | ŷ_A + ŷ_B - μ | r = observed - expected | 표준 |
| Multiplicative | ŷ_A × ŷ_B / μ | r = observed/expected - 1 | 로그 스케일에서 가법과 동치 |
| Product neutrality | fitness_A × fitness_B | r = fitness_AB/(fitness_A·fitness_B) - 1 | Valenzuela 2025 |

**에피스타시스 분류 (5-class, GEARS 기준):**
1. Synergy: r > 0 (관찰 > 기대, 같은 방향 강화)
2. Suppression: r < 0 (관찰 < 기대, 한쪽이 다른 쪽 억제)
3. Neomorphism: 방향이 다름 (새로운 표현형 발현)
4. Buffering: r ≈ 0 (조합 효과 ≈ 더 큰 단일 효과)
5. Additive: |r| < threshold (유의하지 않은 잔차)

**탐지 알고리즘:**
1. 3공식 각각으로 에피스타시스 점수 계산
2. 공식 간 불일치 측정 → "formula epistemic uncertainty"
3. 3공식 모두 유의한 상호작용 → "high-confidence epistasis"
4. ICM 위반 점수로 가중: 높은 ICM 위반 = 높은 에피스타시스 신뢰도

**평가:**

| 지표 | Trivial baseline | 타겟 | 평가 방법 |
|------|-----------------|------|----------|
| Epistasis AUROC | ~0.5-0.6 | > 0.75 | Norman 104 double KO, GEARS GI 라벨 |
| Epistasis Precision | ~0.3 | > 0.6 | Top-k precision at various k |
| GI subtype F1 (5-class) | — | > 0.5 | GEARS 5-subtype 분류와 비교 |
| Formula robustness | — | 3공식 교차 일치율 > 60% | 공식 간 탐지 일치도 |

### 4.6 Phase 6: 능동학습 시뮬레이션 (RQ4)

**설계:**
- 초기: Norman single KO만 관측 가능하다고 가정
- 순차적으로 double KO를 선택하여 관측 (총 104개 중 선택)
- 비교: 불확실성 기반 선택 vs 랜덤 선택 vs 엔트로피 기반 선택

**획득 함수:**
```
α(comb) = U(comb) + λ · diversity(comb, S_selected)
```
- U(comb): Phase 4의 불확실성 점수
- diversity: 이미 선택된 조합과의 다양성 (z_tx 공간 거리)
- λ: 탐색-활용 균형 파라미터

**평가:**
- Top-k recall: 처음 k개 선택에서 에피스타시스 탐지율
  - 랜덤: k/104 (baseline)
  - 타겟: > 2× random
- 학습 곡선: 관측 수 대비 예측 R2

---

## 5. 베이스라인 비교

### 5.1 에피스타시스 탐지 (RQ1, RQ3)

| 계층 | 방법 | 구현 | 예상 성능 |
|------|------|------|----------|
| Trivial | 가법 잔차 = 에피스타시스 | |r_add| > threshold | Recall ~100%, Precision ~0.3 |
| Baseline 1 | Pacalin 방식 | 가법 합 + p<0.05 | 중간 정밀도 |
| Baseline 2 | GEARS | 학습 데이터 GI 분류 | 40% higher precision |
| 소거 | FCR (no ICM) | ICM 정규화 제거 | 불확실성 분해 불가 |
| **제안** | FCR-ICM 잔차 분해 | 3-공식 + ICM 가중 | — |

### 5.2 불확실성 정량화 (RQ2)

| 계층 | 방법 | 구현 | 예상 성능 |
|------|------|------|----------|
| Trivial | MC Dropout 분산 | 30회 forward pass | rho ~0.3-0.5, 미보정 |
| Baseline 1 | Deep Ensemble | 5 모델 앙상블 | rho ~0.4-0.6 |
| **제안** | ICM + MC + 잔차 결합 | 다중 신호원 + 보정 | rho > 0.6 |

### 5.3 능동학습 (RQ4)

| 계층 | 방법 | 구현 | 예상 성능 |
|------|------|------|----------|
| Trivial | 랜덤 조합 선택 | 무작위 샘플링 | k/104 recall |
| Baseline 1 | 불확실도 순 | U(comb)만 사용 | 중간 |
| **제안** | ICM + 잔차 기반 | U + diversity | > 2× random |

---

## 6. 소거 실험 매트릭스

| ID | FCR | ICM | 3-공식 | 잔차 분해 | UQ 결합 | AL | 예상 결과 |
|----|-----|-----|--------|----------|---------|-----|----------|
| A1 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Full model |
| A2 | ✓ | — | ✓ | ✓ | ✓ | ✓ | ICM 정규화 기여 측정 |
| A3 | ✓ | ✓ | — (단일 공식) | ✓ | ✓ | ✓ | 공식 민감도 기여 측정 |
| A4 | ✓ | ✓ | ✓ | — (trivial 분해) | ✓ | ✓ | 잔차 분해 기여 측정 |
| A5 | ✓ | ✓ | ✓ | ✓ | — (MC only) | ✓ | 다중 신호원 기여 측정 |
| A6 | ✓ | ✓ | ✓ | ✓ | ✓ | — | AL 기여 측정 |
| B1 | CPA | — | — | — | MC only | — | CPA baseline |
| B2 | GEARS | — | — | — | — | — | GEARS baseline |

---

## 7. 데이터 확보 및 전처리

| 데이터셋 | 용도 | 접근성 | 전처리 |
|----------|------|--------|--------|
| Norman et al. (2019) | RQ1-4 주 검증 | GEO GSE133344 | Scanpy 표준, 2000 HVG |
| Replogle et al. (2022) | RQ2 교차 세포 유형 UQ | GEO GSE142398 | K562 + RPE1 분리 |
| PORTAL (Tang & Norman, 2026) | RQ3 대규모 검증 (선택) | 접근성 확인 필요 | — |
| MSigDB Hallmark | 경로 평가 기준 | 공개 | 50 gene sets |

**Norman 전처리 파이프라인:**
1. Raw count → Scanpy 일관된 파이프라인
2. 미토콘드리아/리보좀 유전자 필터
3. 2000 HVG 선택
4. 단일 섭동 131개 + 이중 섭동 104개 분리
5. Control 셀 분리 (기준 분포)
6. 정규화 + 로그 변환

---

## 8. 평가 지표 요약

| RQ | 지표 | Baseline | 타겟 | 데이터 |
|----|------|----------|------|--------|
| RQ1 | 잔차 분해 정확도 (합성) | — | r > 0.7 | 합성 + Norman |
| RQ2 | U-Error Spearman rho | ~0.3-0.5 | > 0.6 | Norman + Replogle |
| RQ2 | Coverage (90% CI) | — | 0.85-0.95 | Norman + Replogle |
| RQ3 | Epistasis AUROC | ~0.5-0.6 | > 0.75 | Norman |
| RQ3 | Epistasis Precision | ~0.3 | > 0.6 | Norman |
| RQ3 | GI subtype F1 (5-class) | — | > 0.5 | Norman |
| RQ4 | Top-k recall (에피스타시스) | k/n (random) | > 2× random | Norman AL 시뮬레이션 |

---

## 9. 논문 구조 제안

1. **Introduction**: Perturb-seq 조합 예측의 근본적 한계 — "예측이 실패할 때 무엇을 배울 수 있는가?"
2. **Related Work**: 섭동 예측(GEARS, CPA), 에피스타시스 정의(Valenzuela, Chitra, Ajmal), 인과 표현 학습(FCR, IEM), 불확실성 정량화(CIPHER)
3. **Methods**: FCR-ICM 잔차 분해 → 불확실성 정량화 → 3-공식 에피스타시스 탐지 → 능동학습
4. **Results**:
   - 4.1 RQ1: 잔차 분해가 에피스타시스 신호를 정밀하게 분리하는가?
   - 4.2 RQ2: 다중 신호원 불확실성이 예측 오차와 상관하는가?
   - 4.3 RQ3: 3-공식 에피스타시스 탐지의 정밀도와 견고성
   - 4.4 RQ4: 불확실성 기반 능동학습의 실험 효율 개선
5. **Discussion**: 잔차를 신호로 해석하는 패러다임 전환, 한계(식별 가능성, 공식 민감도), 확장(PORTAL 대규모 검증)

---

## 10. 핵심 리스크와 완화

| 리스크 | 확률 | 영향 | 완화 |
|--------|------|------|------|
| 잔차 분해의 식별 가능성 부족 | 중간 | 높음 | 합성 데이터로 검증; 약화된 식별 조건(IEM) 사용 |
| 에피스타시스 AUROC < 0.75 | 중간 | 높음 | 공식 민감도를 결과로 보고; "high-confidence epistasis" 하위 집합으로 타겟 완화 |
| ICM 위반이 에피스타시스와 무관 | 낮음-중간 | 높음 | 통제 변수 분석; Norman에서 ICM 위반 패턴 검사 |
| MC Dropout + ICM 결합이 단일 신호보다 나쁨 | 낮음 | 중간 | 소거 실험 A5로 검증; 단일 신호로 폴백 |
| Norman 규모 제한 (104 double KO) | 높음 | 중간 | PORTAL로 확장 검증; 교차 검증 활용 |
| 3-공식 결과가 완전히 불일치 | 낮음 | 중간 | 불일치 자체를 "formula epistemic uncertainty"로 보고; Valenzuela 근거로 Product 우선 |

---

## 11. 이전 run 실패에서의 설계 원칙

| 실패 | Run | 교훈 | 본 설계 반영 |
|------|-----|------|-------------|
| R2 기반 평가로 novelty 불충분 | run_02/08 | R2 개선만으로는 publishable하지 않음 | 평가 척도를 Precision/Recall/AUROC로 전환 |
| Mean Shift baseline이 강함 | run_08 | 강한 trivial baseline 존재 시 차별화 어려움 | "잔차=에피스타시스" trivial baseline은 precision 낮음 (0.3) |
| 경쟁 밀도 5개 | run_05 | 포화 시장에서 novelty 불가 | 직접 경쟁자 1개(CIPHER)인 에피스타시스 분해로 전환 |
| Comp loss 불필요 | run_06 | Norman 실데이터에서 comp loss가 도움 안 됨 | 조합 예측은 가법 기준, 복잡한 loss 불필요 |
| 잠재공간-유전자공간 갭 | run_05 | 잠재공간 평가만으로 불충분 | 모든 평가를 유전자 공간에서 수행 |

---

## 12. 컴퓨팅 자원 추정

| 단계 | 예상 시간 | 비고 |
|------|----------|------|
| Norman 전처리 | 1-2h | Scanpy, 메모리 16GB+ |
| FCR-ICM 학습 | 2-4h/실험 | GPU 필수, RTX 4060 Ti |
| 잔차 분해 + UQ | 1-2h | MC Dropout 30회, inference |
| 3-공식 에피스타시스 탐지 | 1h | 공식별 계산 |
| AL 시뮬레이션 | 2-4h | 순차 선택 반복 |
| 소거 실험 (8 구성) | 1-2일 | × 3-fold CV |
| Replogle 교차 검증 (RQ2) | 4-6h | 추가 데이터 처리 |

**총 예상: 3-5일 (GPU 포함)**
