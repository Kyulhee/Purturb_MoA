# FCR-ICM Novelty Assessment — 논문 출판 가능성 평가

**Date:** 2026-04-29 | **Context:** run_09 완료 후, 현재 결과로 논문이 가능한지 평가

---

## 1. 경쟁 논문 지도 (2025-2026 신규 포함)

### 직접 경쟁자 (잔차→에피스태시스 분해)

| 논문 | 방법 | 잔차→에피스태시스 | ICM | UQ | AL | 교세포 전이 |
|------|------|:---:|:---:|:---:|:---:|:---:|
| **FCR-ICM (우리)** | VAE z_x/z_t/z_tx + ICM MMD | ✅ | ✅ | ⚠️ (0.401) | ✅ (5.0x) | ✅ |
| GEARS (Roohani 2023) | GNN + 지식그래프 | ❌ | ❌ | ❌ | ❌ | ❌ |
| CPA (Lotfollahi 2023) | 가법 임베딩 | ❌ | ❌ | ❌ | ❌ | ❌ |
| Pacalin 2025 (Nature Biotech) | CRISPRai 이중섭동 | PARTIAL* | ❌ | ❌ | ❌ | ❌ |
| PORTAL 2026 | Reporter pooled genetics | ❌ | ❌ | ❌ | ❌ | ❌ |
| **C3TL** (arXiv:2603.13051) | 인과 불변성 전이 | ❌ | ✅* | ❌ | ❌ | ✅ (bulk) |

*Pacalin은 단순 가법 합 잔차 사용. 훈련된 예측기 잔차 아님.
*C3TL은 인과 불변성 원리로 cross-domain 전이. 단 bulk 데이터, single-cell VAE 분해 아님.

### 간접 경쟁자 (분해 표현 학습)

| 논문 | 방법 | 표현 분해 | 도메인 불변 | 교세포 전이 | 에피스태시스 | UQ |
|------|------|:---:|:---:|:---:|:---:|:---:|
| **scDRP** (Sun 2025, bioRxiv) | β-VAE + sparsity | ✅ (2-way) | ❌ | ❌ | ❌ | ❌ |
| **SCCVAE** (Liu 2026, PLOS CB) | 인과모델 + VAE | PARTIAL | ❌ | ❌ | ❌ | ❌ |
| **BuDDI** (Davidson 2025, PLOS CB) | VAE 4-분해 | ✅ (4-way) | ✅ (domain inv.) | PARTIAL | ❌ | ❌ |
| **scREPA** (Wang 2026) | Cycle-consistent | ❌ | ❌ | ❌ | ❌ | ❌ |
| **C3TL** (arXiv:2603.13051) | 인과 불변성 | ❌ | ✅ (causal inv.) | ❌ | ❌ | ✅ (bulk) |

**핵심: 잔차→에피스태시스 분해를 하는 논문 = 0개.** 이것이 우리의 가장 큰 novelty. C3TL이 ICM과 유사한 인과 불변성 원리를 사용하나, bulk 데이터에 적용하고 VAE 분해/에피스태시스 없음.

---

## 2. FCR-ICM의 독자적 기여 (Novelty 강도별)

### 강한 novelty (경쟁자 0-1명)

1. **예측 잔차의 체계적 에피스태시스 분해** (RQ1)
   - r = r_epistasis + r_model + r_noise
   - 합성 데이터에서 완벽 검증 (r=1.000)
   - 직접 경쟁자 0개. Pacalin(2025)이 최근접이나 훈련된 예측기 아님
   - **Novelty: ★★★★★**

2. **ICM 정규화로 교세포 zero-shot 전이** (RQ3 from earlier runs)
   - MMD로 z_tx를 cell type 불변으로 만들어 전이 R2 -0.30→0.92
   - BuDDI가 domain invariant 분해를 하나, 목적(bulk→cell-type)이 다르고 불변성 메커니즘이 다름
   - **C3TL (arXiv:2603.13051)이 독립적으로 인과 불변성 원리로 전이를 수행** — 직접 선행. 단 C3TL은 bulk 데이터, single-cell VAE 분해 없음
   - **Novelty: ★★★☆☆** (C3TL이 인과 불변성 전이라는 아이디어 자체는 선행. 우리의 차별화: single-cell VAE + ICM 정규화 메커니즘 + z_tx 분해 구조)

3. **3-공식 에피스태시스 민감도 분석** (RQ3)
   - 가법/곱법/Product neutrality 공식 비교, 방향 일치도 76.6%
   - Valenzuela(2025)가 product neutrality 이론 제시했으나, 예측 잔차에 적용한 것은 우리가 최초
   - Chitra(2025)가 공식 선택의 영향 입증했으나, 섭동 예측 잔차 컨텍스트 아님
   - **Novelty: ★★★★☆**

4. **OOD distance 기반 불확실성 신호** (RQ2)
   - z_tx 공간에서 조합의 "새로움"이 예측 오차와 상관 (rho=0.385)
   - CPA가 임베딩 거리를 불확실성 프록시로 사용한 것이 선행
   - 단, 우리는 ICM 정규화 후의 z_tx 공간에서 측정 → ICM과의 시너지가 새로움
   - **Novelty: ★★★☆☆** (CPA가 선행, ICM+OOD 결합은 새로움)

### 약한 novelty (경쟁자 다수)

5. **능동학습으로 조합 실험 우선순위 결정** (RQ4)
   - 5.0x 개선율 (에피스태시스 점수), 2.6x (OOD)
   - AL for experiment design 자체는 새로운 개념이 아님
   - 단, 섭동 생물학에서 조합 실험 AL은 거의 없음
   - **Novelty: ★★★☆☆**

6. **VAE 표현 분해** (z_x/z_t/z_tx)
   - CPA, BuDDI, scDRP 모두 유사 분해 수행
   - FCR의 3-way 분해 자체는 독창적이나, 분해 아이디어는 일반적
   - **Novelty: ★★☆☆☆**

---

## 3. 약점 분석

### 치명적 약점

| 약점 | 심각도 | 이유 | 완화 가능? |
|------|--------|------|-----------|
| RQ2 rho=0.401 < 0.6 | **HIGH** | UQ가 핵심 기여라면 미달. 단, UQ가 주 기여가 아닐 수 있음 | ⚠️ Replogle에서 재평가 |
| RQ3 AUROC=1.0 순환 | **MEDIUM** | 가법 잔차→effect_size→가법 잔차. Product-GT rho=0.620으로 대체 가능 | ✅ OOD/ICM 점수로 교체 |
| Norman 단일 세포유형 | **MEDIUM** | ICM 진가(다세포유형)를 run_07에서만 부분 검증 | ⚠️ Replogle 이중 세포유형 |
| 에피스태시스 분해 실데이터 검증 불가 | **HIGH** | Ground truth 성분이 없어 직접 검증 불가 | ⚠️ PORTAL 665K 쌍으로 간접 |

### 비치명적 약점

| 약점 | 심각도 | 대응 |
|------|--------|------|
| 가법 기준 R2=0.968로 잔차 작음 | LOW | 잔차가 작다는 것 자체가 발견 (에피스태시스가 약함) |
| GEARS GI 라벨 불가 | LOW | 통계적 ground truth로 대체 |
| comp loss 불필요 (실제데이터) | LOW | 합성에서는 유효. 실제에서 불필요함도 발견 |

---

## 4. 논문 포지셔닝 옵션

### 옵션 A: "잔차 분해 + ICM 전이" 메인 (추천)

**타이틀 방향:** "Decomposing Prediction Residuals into Epistasis and Model Error via Factorized Causal Representations"

**핵심 내러티브:**
1. FCR 분해 → 교세포 zero-shot 전이 (ICM)
2. 잔차 → 에피스태시스 분해 (합성 검증)
3. 3-공식 민감도 분석 → 에피스태시스는 스펙트럼
4. OOD distance가 불확실성의 가장 강한 신호
5. AL로 조합 실험 5.0x 개선

**강점:** novelty가 가장 높은 "잔차→에피스태시스 분해"를 메인으로
**약점:** RQ2 미달을 한계로 서술. 실데이터 에피스태시스 분해 검증 불가를 인정

### 옵션 B: "ICM 교세포 전이" 메인

**타이틀 방향:** "Zero-Shot Cross-Cell-Type Perturbation Prediction via ICM-Regularized Factorized Representations"

**핵심 내러티브:**
1. FCR + ICM → 교세포 전이 R2 -0.30→0.92 (결정적 기여)
2. ICM이 인코더를 더 선형적으로 만들어 조합성 보존
3. 잔차 분해 + 에피스태시스는 확장 응용

**강점:** 가장 강력한 정량 결과 (R2 극적 개선)
**약점:** BuDDI와 부분 중복. 교세포 전이만으로는 novelty 부족

### 옵션 C: "에피스태시스 스펙트럼 + AL" 메인

**타이틀 방향:** "Epistasis as a Spectrum: From Prediction Residuals to Active Combination Experiment Design"

**핵심 내러티브:**
1. 에피스태시스는 이진이 아니라 스펙트럼 (128/128 double-KO가 가법에서 벗어남)
2. 3-공식 민감도 → 76.6% 일치, 공식 선택이 중요하나 치명적이지 않음
3. OOD distance로 불확실한 조합 식별
4. AL로 5.0x 실험 효율 개선

**강점:** 에피스태시스 스펙트럼 관점이 새로움. AL 실용적 기여
**약점:** RQ2 미달. 순환 평가 문제. 이론적 깊이 부족

---

## 5. 판정: 현재 상태로 논문 가능한가?

### 최소 출판 기준 충족 여부

| 기준 | 상태 | 비고 |
|------|------|------|
| Novelty (독자적 기여) | ✅ | 잔차→에피스태시스 분해: 경쟁자 0 |
| Soundness (방법론 건전성) | ⚠️ | 합성 검증 완료, 실데이터 직접 검증 불가 |
| Significance (중요성) | ✅ | 교세포 전이 + 에피스태시스 + AL |
| Completeness (완결성) | ⚠️ | RQ2 미달, RQ3 순환, 단일 데이터셋 |
| Reproducibility | ✅ | 공개 데이터, 코드 완비 |

### 결론

**옵션 A로 top venue 가능성: 35-45%** (NeurIPS/ICML workshop, ISMB, RECOMB)

- 긍정: novelty가 확실 (잔차→에피스태시스 분해, ICM 전이). 알고리즘/응용 벤치마크 충족
- 부정: RQ2 미달, 실데이터 에피스태시스 검증 불가, 단일 데이터셋, C3TL이 인과 불변성 전이 아이디어 선행
- 부정: RQ2 미달, 실데이터 에피스태시스 검증 불가, Norman 단일 데이터셋

**논문 가능하게 만드는 추가 실험 (우선순위순):**

1. **Replogle 다세포유형에서 에피스태시스 파이프라인 재실행** (가장 중요)
   - ICM 진가 재검증 + 더 큰 잔차에서 RQ2 개선 가능성
   - 예상 소요: ~1 run

2. **RQ3 순환 평가 해소** — OOD/ICM 기반 점수로 에피스태시스 탐지
   - 순환적 AUROC를 비순환적 평가로 교체
   - 예상 소요: run_09c 수정

3. **PORTAL 2026 대규모 검증** (가장 영향력 있으나 데이터 접근 필요)
   - 665,856 쌍으로 precision/recall 평가
   - 예상 소요: 데이터 확보 후 ~1-2 run

---

## 6. 경쟁 대비 차별화 매트릭스

| 기능 | FCR-ICM | GEARS | CPA | scDRP | BuDDI | SCCVAE | C3TL |
|------|---------|-------|-----|-------|-------|--------|------|
| 표현 분해 | 3-way | ❌ | 2-way | 2-way | 4-way | PARTIAL | ❌ |
| ICM 정규화 | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ✅* |
| 교세포 전이 | ✅ zero-shot | ❌ | ❌ | ❌ | PARTIAL | ❌ | ✅ (bulk) |
| 잔차→에피스태시스 | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| 3-공식 민감도 | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| 불확실성 정량화 | ⚠️ weak | ❌ | PARTIAL | ❌ | ❌ | ❌ | ❌ |
| 능동학습 | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| 조합 예측 | ✅ R2=0.88 | ✅ | ⚠️ | ⚠️ | ⚠️ | ✅ | ❌ |

*C3TL은 인과 불변성 원리 사용하나 MMD 정규화 방식 아님

**FCR-ICM만이 6개 기능을 동시에 제공.** 이 "통합성" 자체도 novelty.
