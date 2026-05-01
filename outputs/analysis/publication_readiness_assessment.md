# 출판 가능성 면밀 검토

**Date:** 2026-05-01 | **Context:** BioEval + FCR-ICM 전체 프로젝트 출판 가능성 평가

---

## 0. 핵심 구조: 두 개의 독립적 연구 스레드

현재 프로젝트에는 **서로 다른 두 연구**가 혼재:

| | FCR-ICM (run_01-11) | BioEval (run_12-15) |
|---|---|---|
| 질문 | 섭동 예측을 어떻게 개선할까? | 섭동 예측을 어떻게 평가할까? |
| 핵심 주장 | ICM 전이 + 잔차→에피스태시스 분해 | MSE/R2 순위 반전 + BioEval>downstream 상관 |
| novelty_ledger | 없음 (별도 평가) | objects/current/novelty_ledger.yaml |
| result_card | 포함 안 됨 | objects/current/result_card.yaml |
| 독립 출판 가능? | 별도 검토 필요 | **가능** |

**중요:** BioEval이 현재 result_card/validation_readiness_card에 등록된 활성 프로젝트. 아래 평가는 BioEval을 주대상으로 하되, FCR-ICM도 별도로 평가.

---

## 1. BioEval — 출판 가능성 평가

### 1.1 가설별 출판 준비도

#### H1: MSE/R2와 BioEval 방향 지표 간 모델 순위 불일치 (tau < 0.7)

| 증거 유형 | 결과 | 강도 |
|-----------|------|------|
| 시뮬레이션 (3 데이터셋, 11 모델) | K562 tau=0.382, RPE1 tau=0.600, Norman tau=-0.200 | 강함 |
| 실제 Ridge LOO (Norman, 9 모델) | tau(MSE,Dir_deg)=0.500 (PARTIAL), tau(Pearson,Dir_deg)=-0.167 (REVERSAL) | 중간 |
| 교세포 일관성 | 모든 지표 tau > 0.78 (K562↔RPE1) | 강함 |
| S1 민감도 | DEG 임계값 스윕 PASS | 강함 |

**출판 준비도: 75%**

- **강점:** 이중 증거(시뮬+실제). 교세포 일관성으로 체계적 현상 입증. Mean-effect trap 메커니즘 명확히 규명.
- **약점:** 실제 모델 결과가 Norman 1개 데이터셋. K562/RPE1은 one-hot LOO 퇴화로 부적합. 시뮬레이션 tau(-0.20) vs 실제 tau(0.50)의 차이에 대한 설명 부족.
- **리뷰어 반박 가능성:** "9개 모델로 순위를 매기는 것의 통계적 분산이 크다. N=9에서 tau=0.5가 유의한가?" → **Bootstrap CI가 필수**

#### H2: BioEval이 MSE보다 downstream 생물학적 유용성을 더 잘 예측한다

| 증거 유형 | 결과 | 강도 |
|-----------|------|------|
| 시뮬레이션 (3 데이터셋, 90 테스트) | 80/90 pass (88.9%) | 강함 |
| 실제 Ridge LOO (Norman, 15 테스트) | 15/15 pass (100%) | 강함 |
| 핵심 비교 | rho(MSE,dir_disc)=0.678 vs rho(Dir_deg,dir_disc)=1.000 (diff=+0.322) | 강함 |

**출판 준비도: 80%**

- **강점:** H1과 H2의 결합이 매우 강력. "MSE와 BioEval이 다른 모델을 선택하고(H1), BioEval의 선택이 생물학적 유용성과 더 일치한다(H2)"는 논증은 강력한 내러티브.
- **약점:** downstream 과업의 정의(dir_discovery, f1@50)가 "BioEval 지표 자체의 변형"이라는 순환성 의심 가능. → downstream 과업이 BioEval 지표와 독립적인지 명확히 해야.
- **리뷰어 반박 가능성:** "dir_discovery가 방향 정확도 기반이면 순환 아닌가?" → **downstream 과업의 독립성 정당화 필수**

#### H3: BioEval 하에서 trained > baseline인 섭동 체계 존재

| 증거 유형 | 결과 | 강도 |
|-----------|------|------|
| Norman Ridge | 모든 지표에서 trained > true baselines | 중간 |
| Oracle baseline 문제 | raw에서 trained_wins=False → corrected에서 True | 주의 필요 |

**출판 준비도: 50%**

- **강점:** Ridge가 mean_predictor/mean_effect를 압도 (Dir_deg 0.986 vs 0.000/0.571).
- **약점:** Ridge는 선형 모델. "DL > baseline" 주장은 검증 안 됨. Oracle baseline 분류가 ad-hoc.
- **리뷰어 반박 가능성:** "선형 모델이 baseline을 이기는 것은 당연하지 않은가? 논문의 핵심 주장이 아닌 것 같다." → **H3는 보조 주장으로 격하, 또는 DL 모델 확보 후 강화**

### 1.2 Novelty 평가 (BioEval)

novelty_ledger의 분석에 근거:

| 구성요소 | 독립 신규성 | 선행 예술 | 평가 |
|----------|:----------:|----------|------|
| 유전자×섭동 수준 방향 정확도 (BioEval-Dir) | ✅ | Shesha(섭동 수준만), SCALE PDCorr(유전자 분해 없음) | **강함** |
| 효과크기 보정 분석 (BioEval-Cal) | ✅ | 없음 | **매우 강함** |
| 지표-순위 반전 정량화 (Kendall τ) | ✅ | 없음 | **매우 강함** |
| 지표-downstream 과업 상관 분석 | ✅ | 없음 | **매우 강함** |
| DEG precision/recall 곡선 | ❌ | AUPRC/Zhu(2025) | 약함 |
| 최소 보고 기준 제안 | ❌ | DOME/Walsh(2020) | 약함 |
| 베이스라인 위기 원인 판별 | ❌ | Ahlmann-Eltze(2025), Csendes(2025) | 약함 |

**핵심 novelty: 4개 구성요소가 선행 예술 없이 독립적으로 신규.** 이것은 FCR-ICM(아이디어 독창성 ★★★☆☆)과 대조적으로 매우 유리.

### 1.3 시의성 (Timeliness)

| 요소 | 상태 |
|------|------|
| Ahlmann-Eltze (2025) Nature Methods | "DL≤baseline" 위기 진단 — 해결책 없음 |
| Wei et al. (2026) Nature Methods | 27방법 벤치마크 — 순위 반전 분석 없음 |
| SCALE (2026) | "MSE가 mean-effect trap 유발" 지적 — 정량적 분석 없음 |
| Virtual Cell Challenge (2025) | 평가 표준화 요구 — 지표 설계 아님 |

**평가:** 현재 분야의 핵심 질문("왜 DL이 baseline을 못 이기는가?")에 대한 답을 제공하는 타이밍. **시의성 매우 좋음.**

### 1.4 경쟁 밀도

직접 경쟁자 5개(Wei, Zhu, SCALE, Shesha, Csendes)가 있으나, **모두가 진단에 머물고 해결책(순위 반전 정량화, downstream 과업 상관)은 제공하지 않음.** 이것이 BioEval의 가장 큰 강점.

### 1.5 치명적 약점과 완화 가능성

| # | 약점 | 심각도 | 완화 가능? | 필요 작업 |
|---|------|--------|:----------:|----------|
| C1 | K562/RPE1 one-hot LOO 퇴화 — 2/3 데이터셋 실제 모델 확인 불가 | **HIGH** | ✅ | 유전자 수준 feature Ridge (1 run 예상) |
| C2 | downstream 과업의 순환성 의심 | **MEDIUM** | ✅ | downstream 과업 정의를 BioEval 지표와 독립적으로 정당화 |
| C3 | N=9 모델 → 순위 통계 분산 큼 | **MEDIUM** | ✅ | Bootstrap CI + 더 많은 모델 추가 |
| C4 | Ridge만 → H3 "DL > baseline" 미검증 | **MEDIUM** | ⚠️ | GEARS/CPA 학습 필요 (비용 큼) |
| C5 | Norman logFC 스케일 불일치 | **LOW** | ✅ | 정규화 (간단) |

**C1이 가장 중요.** 유전자 수준 feature Ridge로 K562/RPE1 결과를 확보하면 증거 기반이 3배 확대.

### 1.6 BioEval 출판 가능성 — 벤뉴별

| 벤뉴 | 가능성 | 근거 | 필요 추가 작업 |
|------|--------|------|---------------|
| **Nature Methods** | 15-25% | 시의성 좋으나 2/3 데이터셋 한계 + DL 모델 없음 | C1 해결 + DL 모델 + 대규모 벤치마크 |
| **PLOS Computational Biology** | 50-60% | 방법론 건전, novelty 확실, 실험 충분 | C1 해결 + Bootstrap CI |
| **Bioinformatics** | 45-55% | 지표 설계 논문으로 적합 | C1 해결 + Bootstrap CI |
| **ISMB/RECOMB** | 40-50% | 컴퓨테이셔널 바이올로지 적합 | C1 해결 + downstream 과업 독립성 정당화 |
| **NeurIPS/ICML workshop** | 30-40% | ML 벤뉴에서는 생물학적 동기가 약함 | DL 모델 결과 필수 |
| **NeurIPS/ICML main** | 10-15% | ML novelty 부족 — "지표 설계"는 생물학 벤뉴에 더 적합 | 근본적 방향 전환 필요 |

**현실적 판단: PLOS CB / Bioinformatics급이 가장 적합. C1만 해결하면 50% 이상 가능.**

---

## 2. FCR-ICM — 출판 가능성 평가 (참고용)

### 2.1 핵심 결과 재확인

| 가설 | 결과 | 강도 |
|------|------|------|
| ICM 교세포 zero-shot 전이 (RQ3) | R2 -0.30→0.92 | 매우 강함 |
| 잔차→에피스태시스 분해 (RQ1) | 합성 r=1.000 | 강함 (합성만) |
| MC Dropout UQ (RQ2) | Replogle rho=0.660 PASS / Norman 0.401 미달 | 중간 |
| 에피스태시스 AL (RQ4) | 5.0x random | 강함 |
| 3-공식 민감도 | 76.6% 일치 | 중간 |

### 2.2 FCR-ICM의 구조적 한계 (structural_failure_analysis에서 확인)

1. **아이디어 독창성 한계:** "잔차→비가법 신호"는 Tan(2018)/Diamond(2024)에서 이미 일반화. MoCHI/Otwinowski가 생물학에서 유사 접근.
2. **실데이터 에피스태시스 분해 검증 불가:** ground truth 성분이 없음. 이것은 **이 연구 방향의 근본적 한계**.
3. **OOD 비견고:** Norman +0.385 vs Replogle -0.340.
4. **단일 아키텍처, 하이퍼파라미터 튜닝 없음.**
5. **C3TL이 인과 불변성 전이 아이디어 선행.**

### 2.3 FCR-ICM 출판 가능성

이전 평가(30-40%)를 유지하되, **하향 요인이 더 많음:**

- novelty_ledger에서 FCR-ICM의 novelty는 ★★★☆☆~★★★★☆로 평가되었으나, 이는 FCR-ICM 스레드의 것이고, BioEval novelty_ledger는 FCR-ICM을 평가하지 않음
- structural_failure_analysis에서 식별된 5개 구조적 문제가 모두 유효
- **가장 큰 문제:** 실데이터 에피스태시스 분해 검증이 불가능하다는 것이 **근본적 한계** (추가 실험으로 해결 불가)

**FCR-ICM 단독 출판 가능성: 25-35%**

---

## 3. 두 스레드의 결합 — 시너지 가능성?

### 3.1 결합 논문 아이디어

> "BioEval: 생물학적 충실도 기반 섭동 예측 평가 지표 — FCR-ICM을 사례 연구로"

- H1/H2를 일반적 프레임워크로 제시
- FCR-ICM을 BioEval로 평가한 결과를 사례로 포함
- 문제: FCR-ICM의 결과가 BioEval의 주장을 돕지 않음 (FCR-ICM은 평가 대상이지 평가 방법이 아님)

### 3.2 결합의 문제

- 두 스레드는 질문이 다름: "어떻게 예측할까" vs "어떻게 평가할까"
- 결합하면 논문 초점이 흐려짐
- **권장: 분리 출판**

---

## 4. 종합 판정

### 4.1 BioEval 출판 — 최종 평가

| 항목 | 평가 | 근거 |
|------|------|------|
| **Novelty** | ★★★★☆ | 4개 구성요소 선행 예술 없음. "지표-순위 반전 정량화"와 "지표-downstream 과업 상관"은 완전 신규 |
| **Soundness** | ★★★★☆ | 합성+실제 이중 검증. 교세포 일관성. S1 민감도. 단 1/3 데이터셋만 실제 모델 확인 |
| **Significance** | ★★★★★ | 분야 핵심 문제(DL≤baseline 위기)에 직결. 시의성 최고 |
| **Completeness** | ★★★☆☆ | C1(2/3 데이터셋 한계), C4(DL 모델 없음)가 미해결 |
| **Reproducibility** | ★★★★★ | 공개 데이터(Norman, Replogle), 코드 완비 |

**종합 출판 가능성: 45-55%** (PLOS CB / Bioinformatics급)

- C1 해결(유전자 수준 feature Ridge) 시 → **55-65%**
- C1 + C4 해결(DL 모델 추가) 시 → **60-70%**

### 4.2 최소 출판 기준 달성 여부

| 기준 | 달성? | 비고 |
|------|:-----:|------|
| Novelty (독자적 기여 ≥1개) | ✅ | 4개 구성요소 완전 신규 |
| Soundness (방법론 건전성) | ✅ | 이중 검증, 민감도 분석 |
| Significance (시의성+중요성) | ✅ | Nature Methods 2025 위기 직결 |
| Completeness (3개 이상 데이터셋) | ⚠️ | 시뮬레이션 3개 OK, 실제 1개만 |
| Reproducibility | ✅ | 완비 |
| 통계적 견고성 | ❌ | Bootstrap CI 없음. N=9 순위 분산 큼 |

**최소 기준: 5/6 달성. Bootstrap CI 추가 시 6/6.**

### 4.3 FCR-ICM 출판 — 최종 평가

**종합 출판 가능성: 25-35%**

- 근본적 한계(실데이터 에피스태시스 검증 불가)가 추가 실험으로 해결 불가
- 아이디어 독창성이 Tan(2018)/Diamond(2024)로 인해 약화
- 교세포 전이 결과(R2 -0.30→0.92)는 강력하나 C3TL 선행
- **권장:** BioEval 출판 후, 보완 작업(FDR 제어, PORTAL 검증)으로 후속 논문 검토

---

## 5. 권장 로드맵

### Phase A: BioEval 최소 출판 (2-3 run 예상)

| 단계 | 작업 | 목적 | 예상 소요 |
|------|------|------|----------|
| A1 | 유전자 수준 feature Ridge (K562/RPE1) | C1 해결 — 3/3 데이터셋 실제 모델 확보 | 1 run |
| A2 | Bootstrap CI for all tau/rho | 통계적 견고성 확보 | 0.5 run |
| A3 | Norman logFC 스케일 보정 | C5 해결 | 0.3 run |
| A4 | downstream 과업 독립성 정당화 | C2 해결 — 순환성 의심 제거 | 문서 작업 |
| A5 | 논문 초안 (Stage 05) | — | — |

**Phase A 완료 후 출판 가능성: 55-65% (PLOS CB / Bioinformatics)**

### Phase B: BioEval 강화 (선택적, 2-3 run)

| 단계 | 작업 | 목적 | 예상 소요 |
|------|------|------|----------|
| B1 | GEARS/CPA DL 모델 학습 | C4 해결 — H3 "DL > baseline" 확정 | 2-3 run |
| B2 | 추가 downstream 과업 (pathway enrichment, gene set overlap) | H2 강화 | 1 run |
| B3 | 더 많은 모델 다양성 (NN, transformer) | 순위 분산 감소 | 1-2 run |

**Phase B 추가 후 출판 가능성: 65-75% (PLOS CB) / 30-40% (Nature Methods)**

### Phase C: FCR-ICM 후속 (별도)

- BioEval 출판 이후 검토
- FDR 제어(knockoff/순열 검정) 추가가 가장 시급
- PORTAL 665K 쌍 검증이 가장 영향력 있으나 데이터 접근성 불확실

---

## 6. 정직한 평가: 위험 요소

### 6.1 BioEval의 숨겨진 위험

1. **downstream 과업 순환성:** dir_discovery와 f1@50이 방향 정확도 기반 → "방향 지표가 방향 기반 과업을 잘 예측한다"는 순환 가능성. 이것이 리뷰에서 가장 공격받을 포인트.
   - **완화:** pathway enrichment recovery 등 방향과 무관한 downstream 과업 추가 필요

2. **Norman 특이성:** Norman의 극적 반전(tau=-0.20)이 낮은 DEG 비율(1.53%)과 적은 섭동 수(283)의 산물일 수 있음. "일반적 현상"인지 "Norman의 특이성"인지 불확실.
   - **완화:** K562/RPE1 실제 모델 결과가 필수 (C1)

3. **Ridge가 "진짜" 모델인가:** Ridge LOO를 "실제 학습 모델"로 분류했으나, 리뷰어가 "이것은 여전히 베이스라인이다"라고 주장 가능.
   - **완화:** GEARS/CPA 등 "진짜 DL 모델" 결과가 최소 1개 필요 (H3용은 아니더라도 H1/H2 확인용)

4. **순위 반전의 실제 영향:** "순위가 다르다"는 것과 "순위가 다르기 때문에 실제로 더 나은 모델을 놓친다"는 것은 다름. 현재 후자의 증거가 H2뿐인데, H2의 downstream 과업이 순환적.
   - **완화:** 실제 실험 설계 시나리오에서 "MSE가 선택한 모델 vs BioEval이 선택한 모델"의 성과 차이를 보여주는 것이 가장 강력

### 6.2 FCR-ICM의 근본적 딜레마

- 기술적 성공(RQ1/RQ3 통과)과 novelty 성공이 분리됨
- "성공의 역설": 성공한 결과(교세포 전이)가 다른 선행(BuDDI, C3TL)과 중복
- 가장 독창적인 기여(잔차→에피스태시스 분해)는 실데이터 검증 불가
- **이 딜레마는 구조적이며, 추가 실험으로 해결 불가**

---

## 7. 결론

### BioEval: 출판 가능, Phase A 완료 후 제출 권장

- 현재 상태로도 **최소 출판 기준 충족** (단 통계적 견고성 Bootstrap CI 필요)
- Phase A (C1 해결 + Bootstrap CI) 완료 시 **PLOS CB / Bioinformatics급 출판 55-65% 가능**
- H1+H2 결합 내러티브가 강력. H3은 보조 주장으로 격하.
- **가장 시급한 작업:** 유전자 수준 feature Ridge (C1 해결) + Bootstrap CI

### FCR-ICM: 단독 출판 어려움, BioEval 이후 재검토

- 근본적 한계(실데이터 검증 불가)가 추가 실험으로 해결 불가
- 25-35% 가능성. FDR 제어 + PORTAL 검증이 전제조건
- BioEval 논문에서 FCR-ICM을 "사례 연구"로 활용하는 것은 가능하나, 결합 출판은 비권장
