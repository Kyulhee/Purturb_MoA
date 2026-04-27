# Stage 02 — Framing

## Loopback 기록
- **2026-04-28**: Analysis → Literature Review 루프백. run_08 기준선 비교에서 Mean Shift R2=0.82, FCR+ICM R2=0.92 (+0.10만 우위)로 교세포 전이 novelty 불충분 판정. 문헌에서 누락된 직접 경쟁자 5개(BuDDI, C3TL, scDRP, XTransferCDR, CPA) 발견. FCR-ICM의 핵심 아이디어(VEA 분해+도메인 불변성+교세포 전이)가 이미 포화 상태. docs/01, docs/02에 신규성 세이프가드 추가 완료.
- **2026-04-27**: Analysis → Framing 루프백. 기존 NAP 프레임워크의 novelty 부족으로 연구 질문 전면 재설계. 대사 네트워크 surrogate → Perturb-seq cross-cell-type combinatorial prediction으로 도메인 전환. 가설 선택 과정은 `outputs/analysis/run_04/cross_domain_novelty_scan.md` 참조.
- **2026-04-27**: Analysis → Framing 루프백. 기존 NAP 프레임워크의 novelty 부족으로 연구 질문 전면 재설계. 대사 네트워크 surrogate → Perturb-seq cross-cell-type combinatorial prediction으로 도메인 전환. 가설 선택 과정은 `outputs/analysis/run_04/cross_domain_novelty_scan.md` 참조.

## Workflow
1. `docs/02_framing.md` 가이드 확인
2. `stages/01_literature_review.md`에서 인사이트 확인
3. 연구 질문, 베이스라인, 타겟 성능 수치 정의
4. 산출물 → `outputs/framing/run_XX/`에 저장
5. 아래 지식 업데이트 (검증된 인사이트만 통합)
6. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 연구 질문이 단일 문장으로 기술 가능
- 베이스라인 수치가 문헌 근거와 함께 명시됨
- 타겟 성능 수치가 베이스라인 기반으로 설정됨
- 평가 지표 및 데이터셋 확정

---

## 가설 선택 과정

### 탐색 범위
Perturb-seq 5대 미해결 문제(outputs/analysis/run_04/perturb_seq_unsolved_problems_review.md) × 10개 교차 도메인 아이디어(outputs/analysis/run_04/cross_domain_novelty_scan.md)에서 3개 고참신성 가설 도출.

### 후보 가설과 선택/기각 근거

| 가설 | 해결 격차 | 참신성 | 실행가능성 | 영향 | 판정 |
|------|----------|--------|-----------|------|------|
| **H1+H2: 인과 조합** | Problem 2+3 | VERY HIGH | MEDIUM | VERY HIGH | **선택** |
| H1: 인과 불변성 only | Problem 3 | HIGH | MEDIUM | HIGH | 기각 — Problem 2 미해결 |
| H2: 조합 모듈 only | Problem 2 | HIGH | MEDIUM | VERY HIGH | 기각 — 교차 세포 유형 전이 불가 |
| H3: 기계적 해석가능성 감사 | Problem 1+4 | VERY HIGH | MEDIUM | HIGH | 기각 — 진단적 기여만, 모델 개선 안 함 |

### H1+H2 결합 선택의 논리
1. **두 최대 격차 동시 해결**: Problem 2(조합 폭발) + Problem 3(교차 세포 유형 전이) — 단일 가설로 미해결 문제 2개 타격
2. **논리적 시너지**: z_tx가 ICM에 의해 인과적이고 불변 → 불변 모듈 구조가 자연스럽게 조합 지원 → 한 세포 유형에서의 조합이 다른 세포 유형에서도 성립
3. **테스트 가능한 예측 연쇄**: (a) z_tx 불변성 검증 → (b) 모듈-경로 정렬 검증 → (c) 조합 예측 검증 → (d) 실패 = 상위성 탐지
4. **대안 가설보다 명확한 우위**: H1만으로는 조합 불가, H2만으로는 전이 불가, H3은 모델 제안 안 함. 결합만이 세 가지를 모두 충족

### 기각된 가설 참조
- H1/H2/H3 개별 상세: `outputs/analysis/run_04/cross_domain_novelty_scan.md` Top 3 section
- 10개 도메인 전체 탐색 결과: 동일 문서 Domain 1-10

---

## 연구 질문

**"Factorized causal representations를 Independent Causal Mechanism 원리로 정규화했을 때, 섭동 효과가 세포 유형에 걸쳐 불변하는 인과 모듈로 분해되며, 이 모듈들의 조합이 타겟 세포 유형의 학습 데이터 없이 조합 섭동 효과를 예측할 수 있는가?"**

### 하위 질문
1. **RQ1 (불변성)**: ICM 정규화된 z_tx가 세포 유형 간에 불변(invariant)인가?
2. **RQ2 (조합성)**: 경로 수준 모듈이 학습 가능한 상호작용 함수로 조합되어 단일 섭동 데이터만으로 조합 효과를 예측할 수 있는가?
3. **RQ3 (제로샷 전이)**: RQ1+RQ2가 결합될 때, 소스 세포 유형의 데이터만으로 타겟 세포 유형의 조합 섭동 효과를 제로샷으로 예측할 수 있는가?

### 도출 근거

**해결하는 격차 (Perturb-seq 미해결 문제)**:
- **Problem 2 — 조합 폭발**: 단일 섭동 데이터만으로 조합 효과를 예측하는 원칙적 프레임워크 부재. GEARS는 조합 학습 데이터 필수, CPA는 OOD 시 DEG 0.85→0.38 급락, SAMS-VAE는 가법 분해 가정이 상잔 작용에 실패
- **Problem 3 — 교차 세포 유형 전이**: GEARS 명시적으로 "not designed for cross-cell-type transfer". Cell-JEPA는 기저 상태 재구성은 개선하나 효과 크기 추정은 불가. 중간층 분석(arXiv:2604.14838)은 최적 임베딩 층이 작업/맥락 의존적임을 시사

**이론적 기반 (교차 도메인 차용)**:
- **FCR** (Mao et al., arXiv:2410.22472): z_x(공변량), z_t(처리), z_tx(상호작용) 분리 + 블록 식별 가능성 증명 → 분해된 인과 표현의 기반 제공
- **ICM** (Reizinger et al., Schölkopf group, arXiv:2406.14302, ICLR 2025): 메커니즘의 자율성 원리 — 인과 메커니즘은 맥락에 독립적 → z_tx의 교차 세포 유형 불변성에 대한 이론적 근거
- **scBIG** (arXiv:2602.04901): 유전자 프로그램 구조가 모듈 표현에 도움 (6.7% 개선) — 하지만 조합 규칙 없음 → 본 가설이 조합 규칙을 제공

**선행 연구와의 차이 (novelty 근거)**:
- FCR은 분해된 표현을 제공하나 교차 세포 유형 전이나 조합 예측에 적용하지 않음
- ICM은 이론 원리이나 섭동 예측의 정규화 항으로 사용된 적 없음
- MapPFN(arXiv:2601.21092)은 인컨텍스트 학습으로 제로샷 예측을 보이나 명시적 조합 구조나 인과 분해 없음
- C3TL(arXiv:2603.13051)은 인과 유도 편향으로 교차 세포 유형 전이를 시도하나 조합 예측 없음

### 질문이 답하는 것과 답하지 않는 것
- **답하는 것**: 섭동 효과의 인과 모듈 분해가 가능한가? ICM이 교차 세포 유형 불변성을 보장하는가? 불변 모듈 조합이 조합 섭동을 예측하는가? 조합 예측 실패가 상위성(epistasis) 탐지로 기능하는가?
- **답하지 않는 것**: 섭동 예측의 정보이론적 가변성 한계 (Problem 4), 측정 필요성-불일치 해소 (Problem 5) → 후속 연구

### 평가 전략
| 평가 방식 | 지표 | 베이스라인 | 타겟 | 비고 |
|-----------|------|----------|------|------|
| 교차 세포 유형 단일 섭동 | R2 (gene-level DEG) | GEARS intra R2~0.6, inter ~0 | R2 > 0.5 | 소스→타겟 제로샷 |
| 조합 섭동 예측 | R2 (double KO) | Additive baseline R2~0.3 | R2 > 0.5 (독립 경로) | 단일 KO 데이터만 사용 |
| 불변성 검증 | z_tx correlation (cross-cell) | Random baseline < 0.2 | r > 0.7 | ICM 정규화 vs 무정규화 |
| 상위성 탐지 | Epistasis recall | — | Recall > 0.6 | 조합 잔차 → 상위성 분류 |
| 모듈-경로 정렬 | ARI (vs MSigDB) | — | ARI > 0.4 | 학습된 모듈이 알려진 경로와 일치 |

### 베이스라인 계층
1. **Additive baseline**: z_A + z_B (선형 조합) — 하한
2. **CPA** (Lotfollahi et al., 2023): 조합 섭동 오토인코더, OOD 시 DEG 0.85→0.38 급락
3. **GEARS** (Roohani et al., 2023): GNN + GRN, 교차 세포 유형 미지원
4. **FCR (no ICM)**: 분해된 표현이지만 불변성 정규화 없음 — 소거 실험
5. **FCR + ICM (본 제안)**: 인과 불변 모듈 + 조합 함수 — 제안 방법

### 데이터
| 데이터셋 | 세포 유형 | 섭동 유형 | 규모 | 용도 |
|----------|----------|----------|------|------|
| Norman et al. (2019) | K562 | 131 single + 104 double KO | ~100K cells | 조합성 검증 (RQ2) |
| Replogle et al. (2022) | K562 + RPE1 | Genome-scale CRISPRi | ~2.5M cells | 교차 세포 유형 전이 (RQ1, RQ3) |
| Dixit et al. (2016) | K562 | 7 TF KO | ~7K cells | 소규모 검증 |
| sci-Plex (Srivatsan et al., 2020) | 3 cell lines | 188 compounds | ~649K cells | 약물 섭동 확장 (선택) |

### 핵심 이론적 전제와 위험
1. **전제**: z_tx가 ICM에 의해 세포 유형 불변 → **위험**: 식별 가능성 가정이 현실 데이터에서 성립하지 않을 수 있음
2. **전제**: 경로 수준 모듈이 조합 가능 → **위험**: 상위성이 비조합적일 수 있음. 단, 부분적 조합 구조조차 가치 있음
3. **전제**: 소스 세포 유형의 모듈이 타겟에 전이 → **위험**: 세포 유형 특이적 경로 활성도 차이가 모듈 구조를 변형할 수 있음. C3TL(arXiv:2603.13051)은 bulk 데이터로 경쟁적 성능을 보여 가능성 시사

---

## 후속 질문 (본 질문 달성 후)

**"인과 모듈 조합의 실패 모드가 생물학적 상위성의 체계적 탐색 도구로 기능할 수 있는가?"**

- 도출 근거: 조합 잔차 = 예측 불가능한 상호작용 = 후속 실험 우선순위
- 전제조건: RQ2 달성 (조합 예측 R2 > 0.5 for independent pathways)
- 대안 가설 참조: `outputs/analysis/run_04/cross_domain_novelty_scan.md` (H3: 기계적 해석가능성 감사, Domain 4: 등변 신경망 등)

---

## Run 이력 (세부 내용은 outputs/framing/run_XX/ 참조)
- **run_01** (2026-04-25): NAP 기반 프레이밍 완료. GNN 가치 조건 C1-C6 정의, 베이스라인 XGBoost R2=0.91 설정
- **run_02** (2026-04-27): NAP novelty 부족으로 루프백. 인과 조합 가설(H1+H2 결합)로 전면 재설계
