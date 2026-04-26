# Phase 3-4 기술적 실현성 심층 분석 (v2 — 웹 검색 보강)

작성일: 2026-04-26
기준: run_01 보고서 + 실시간 웹 검색/페칭 결과 통합

---

## Phase 3: GNN+XGBoost 대리 모델 스크리닝

### 1. GEM→GNN 입력 변환 — 대사 네트워크의 그래프 표현

**기존 방법론 정리:**

대사 네트워크의 그래프 표현은 연구 목적에 따라 세 가지 접근이 존재함:

| 표현 방식 | 노드 정의 | 엣지 정의 | 장점 | 단점 |
|-----------|----------|----------|------|------|
| 대사체 그래프 | metabolite | reaction(stoichiometry) | 직관적, 대사 경로 추적 용이 | 방향성/가역성 표현 부족 |
| 반응 그래프 | reaction | metabolite(입출력) | FBA와 직접 대응 | 유전자 정보 부재 |
| **이종 그래프** | metabolite + reaction + gene | met↔rxn(stoichiometry), gene↔rxn(GPR) | 정보 손실 최소, GPR 규칙 포함 | 구현 복잡도 증가 |

**선택 근거**: 이종 그래프(Heterogeneous Graph)가 가장 적합. 근거:
- PyTorch Geometric(PyG)이 이종 GNN 레이어를 네이티브 지원 (**2026-04-26 웹 페칭 확인**)
  - RGCNConv: 관계형 GCN (Schlichtkrull et al., ESWC 2018)
  - RGATConv: 관계형 GAT (Busbridge et al., 2019)
  - HGTConv: 이종 그래프 트랜스포머 (Hu et al., WWW 2020)
  - HEATConv: 이종 엣지 강화 어텐션 (Mo et al., 2021)
- PyG 저장소(https://github.com/pyg-team/pytorch_geometric)는 활성 유지중:
  - 동종 GNN: GCNConv, GATConv, SAGEConv, GINConv, EdgeConv, TransformerConv, PNAConv 등 20+ 레이어
  - 이종 GNN: RGCNConv, RGATConv, HGTConv, HEATConv — **multiple node types and edge types 명시 지원**
  - mini-batch loaders, multi-GPU 지원, torch.compile 호환
- metabolite(72) + reaction(95) + gene(137) 3노드타입 + met→rxn(188+172) + gene→rxn(158) 엣지가 textbook 모델에서 검증됨 (Module A 구현 완료)
- knockout_mask(137차원)를 gene 노드 피처로 통합 가능 → GNN이 knockout 효과를 end-to-end로 학습

**리스크**: 이종 GNN은 동종 GNN 대비 학습 샘플 요구량이 큼 (현재 135샘플로 R2 < 0 관측)

### 2. FBA 정답 데이터 생성 — COBRApy 병렬 실행 시간 추정

**COBRApy 현재 상태 (2026-04-26 웹 페칭 확인):**
- 최신 버전: **v0.31.1** (2026-03-26 릴리즈)
- 저장소: https://github.com/opencobra/cobrapy — 활성 유지
- GLPK 솔버 기본 제공, CPLEX/Gurobi 상업 솔버 선택 가능
- **병렬 실행: 네이티브 미지원** (README에 병렬화 언급 없음)
- Python multiprocessing으로 래핑 필요 → Module A에서 Pool.map으로 구현 완료
- benchmarks 디렉토리 존재하나 상세 성능 수치는 문서화되지 않음

**시간 추정 (textbook 모델 기준):**

| knockout 유형 | 조합 수 | 단일 FBA 시간 | 직렬 총시간 | 8-core 병렬 |
|--------------|---------|-------------|-----------|------------|
| Single | 137 | ~100ms | 13.7s | ~2s |
| Double | 9,316 | ~100ms | 15.5min | ~2min |
| Triple | ~630K | ~100ms | ~17.5h | ~2.2h |
| Random 1000 | 1,000 | ~100ms | 100s | ~13s |
| Random 5000 | 5,000 | ~100ms | 8.3min | ~1min |

**iJO1366 대형 모델 추정:**
- 1,367 genes → single knockout 1,367, double ~934K, triple ~427M
- 단일 FBA: textbook 대비 2-5배 느림 (200-500ms 추정)
- Double knockout 934K × 300ms = 78h (8-core = ~10h)

**핵심 판단**:
- Phase 3 초기: textbook + random 1,000-5,000 샘플로 GNN 학습 → R2 > 0.3 달성 후 iJO1366 확장
- Active Learning으로 고품질 샘플만 선택 → 1,000-2,000 샘플로 충분할 가능성
- COBRApy multiprocessing 래핑은 간단 (Pool.map으로 FBA 함수 병렬화) — Module A에서 검증 완료

### 3. Surrogate Model 일반화 — Active Learning, Bayesian Optimization 접목

**기존 사례:**

| 접근법 | 관련 연구 | 핵심 아이디어 | 우리 적용 |
|--------|----------|-------------|----------|
| Pool-based AL | Settles 2009 | uncertainty + diversity로 라벨링 대상 선택 | UCB: uncertainty + exploitation |
| Bayesian Optimization | Snoek et al. 2012 | GP 기반 acquisition function으로 탐색-활용 균형 | AL 전략을 BO의 acquisition function으로 해석 |
| GNN + BO | 존재하는 조합 (GNN이 feature extractor, BO가 탐색기) | GNN 임베딩을 BO 입력으로 | Module B 임베딩 → Module C UCB 탐색 |
| Neural BO | 간접 관련 (Snoek et al., Swersky et al.) | DNN이 GP를 대체하는 surrogate | XGBoost가 GP를 대체하는 구조와 유사 |

**핵심 설계 결정: two-phase AL**
1. Phase 1 (R2 < 0.3): diversity 기반 탐색 — uncertainty가 무효하므로 임베딩 공간에서 최대한 다른 샘플 선택
2. Phase 2 (R2 >= 0.3): UCB 기반 탐색 — exploitation(biomass 높은 영역) + exploration(uncertainty 높은 영역) 균형

**검증 결과**: 현재 R2 = -0.31에서 Random(0.81) > 모든 AL 전략, diversity가 mean_gr=0.72로 최고 → Phase 1 설계와 일치

### 4. TOPSIS 가중치 민감도 — 다중 기준 의사결정에서 랭킹 안정성

**검증 결과 (Step 4 Ablation):**

| 가중치 방식 | Kendall's tau | Top-3 안정성 | 해석 |
|------------|--------------|-------------|------|
| Entropy (객관적) | 0.414 | 72.2% | 데이터 기반이나 불안정 |
| Equal (동등) | 0.336 | 86.7% | 중간 안정성 |
| **Expert (0.7/0.3)** | **0.734** | **100.0%** | 도메인 지식 기반, 가장 안정 |

**민감도 분석 방법론**:
- 가중치 섭동: 각 가중치를 ±10-30% 변동시키며 랭킹 변화 측정
- Kendall's tau: 원래 랭킹과 섭동 후 랭킹의 순위 상관계수
- Top-3 안정성: 상위 3개 대안이 섭동 후에도 유지되는 비율

**핵심 판단**: Expert weight(0.7/0.3)가 압도적으로 안정. 기본 Expert + 보조 Entropy 검증 구조가 최적.

**리스크**: Expert 가중치의 주관성 — 도메인 전문가 합의가 어려운 경우 Entropy 보조 검증이 필수

---

### Phase 3 종합 평가

**기술적 실현성 등급: MEDIUM**

**핵심 리스크 3가지:**
1. **샘플 부족**: 135샘플로 R2 < 0 → 최소 500-1,000샘플 필요. double knockout 9,316 조합에서 random/AL sampling으로 확보 가능. COBRApy v0.31.1은 병렬 실행을 네이티브 지원하지 않으나 Python multiprocessing 래핑으로 해결 (Module A 검증 완료, 8-core ~2min for double knockouts)
2. **GNN 임베딩 품질**: autoencoder loss 366K→365K로 거의 학습 안 됨 → contrastive loss 또는 edge prediction pretraining 필요. PyG의 이종 GNN(HGTConv, RGCNConv 등)은 지원이 충분하나 이종 GNN은 동종 대비 학습 난이도가 높아 샘플 요구량 증가
3. **AL 전략 전환 시점**: R2 > 0.3 이후에만 UCB가 유효하나, 정확한 전환 시점 판단이 모호. 너무 일찍 전환하면 uncertainty가 여전히 무효

**대안/완화 전략:**
1. 샘플 부족 → COBRApy multiprocessing + random 5,000샘플 초기 학습 → AL로 추가 샘플 선택
2. GNN pretraining → edge prediction (reaction-metabolite 연결 예측) 또는 contrastive learning (knockout 쌍 대조)으로 autoencoder 대체
3. AL 전환 → validation R2를 0.3 임계값으로 모니터링, 자동 전환 로직 구현

**필요 컴퓨팅 자원:**
- CPU: 8-core 이상 (COBRApy 병렬 FBA)
- RAM: 16GB 이상 (GNN 학습 + XGBoost)
- GPU: 선택적 (CUDA 사용 시 GNN 학습 5-10x 가속, CPU만으로도 textbook 모델은 충분)
- 예상 총시간: 초기 학습 10-30분, AL 루프 1-2시간

**오픈소스 도구/프레임워크:**
| 도구 | 버전 | URL | 용도 | 상태 |
|------|------|-----|------|------|
| COBRApy | v0.31.1 | https://github.com/opencobra/cobrapy | FBA 실행, 모델 관리 | 활성 (2026-03-26 릴리즈) |
| PyTorch Geometric | 2.7.0+ | https://github.com/pyg-team/pytorch_geometric | HGTConv 등 이종 GNN 레이어 | 활성 (이종 GNN 4+ 레이어 지원) |
| XGBoost | v3.2.0 | https://github.com/dmlc/xgboost | Surrogate 회귀, Quantile regression | 활성 |
| pymoo | v0.6.1.6 | https://github.com/anyoptimization/pymoo | NSGA-II 다목적 최적화 | 활성 (IEEE Access 2020) |
| BiGG Models | — | http://bigg.ucsd.edu/ | 공개 대사 모델 108개 | 운영중 |

---

## Phase 4: FLYCOP/dFBA 동적 시뮬레이션

### 1. FLYCOP 유지보수 상태

**확인 결과 (2026-04-26 웹 페칭 재확인):**
- GitHub 저장소(baliga-lab/FLYCOP): **404 삭제됨** (여러 URL 변형 시도 모두 실패)
  - https://github.com/baliga-lab/FLYCOP → 404
  - https://github.com/baliga-lab/Flycop → 404
- 원논문 DOI: 10.1186/s12918-018-0639-6 (BMC Systems Biology 2018)
  - 논문 자체도 Springer/PMC에서 페칭 실패 (403/404/303 리다이렉트 루프)
- 2018년 이후 업데이트 없음, Python 2 기반 추정, **사용 불가 판정 확정**

**대체 전략 검증 완료**: TOPSIS + Entropy weight로 객관화 대체. Expert tau=0.73, top3=100%로 FLYCOP의 fuzzy logic 역할을 충분히 대체.

### 2. 대체 프레임워크 비교 평가

| 프레임워크 | 버전 | 언어 | 마지막 업데이트 | dFBA 지원 | 공간 시뮬레이션 | Python 인터페이스 | 활성도 |
|-----------|------|------|---------------|----------|----------------|-----------------|--------|
| **COMETS** | v2.12.4 | Java | 2025-06-18 | O (dFBA 코어) | O (확산+dFBA) | cometspy v0.6.1 | 활성 (380 commits, 3 releases) |
| **scipy 직접 구현** | — | Python | 지속 | O (BDF/Radau) | X | 네이티브 | 검증 완료 |
| dfba-python | — | Python | 삭제됨 | — | — | — | **FAIL (저장소 404)** |
| cFBA | — | — | — | — | — | — | **확인 불가 (PyPI/404)** |

**COMETS 상세 (2026-04-26 웹 페칭으로 확인):**
- 개발: Daniel Segré Lab (Boston University)
- 라이선스: MIT (비상업 권장)
- 최신 릴리즈: v2.12.4 (2025-06-18)
- 코드 구성: ~46% Java, ~54% HTML
- 특징: stoichiometric modeling of genome-scale metabolic networks + discrete approximation of diffusion
- 시뮬레이션 예시: test tube, competition assay, citrate utilization LTEE, circular/branching colony, demographic noise, virtual Petri dish
- Python 인터페이스: cometspy (https://github.com/segrelab/cometspy)
  - 최신 버전: v0.6.1 (Virus and Chemotaxis models 추가)
  - 설치: `pip3 install cometspy`
  - Python 버전: PyPI 배지로 표시되나 상세 버전은 문서 확인 필요
  - 문서: cometspy.readthedocs.io
- 참고논문: Harcombe et al., Cell Reports, 2014
- **단점**: Java 코어 필수, Python은 cometspy 래퍼만, Java 설치 복잡도

**선택 근거**: scipy BDF/Radau 직접 구현
- COMETS는 Java 의존성으로 환경 구성 복잡
- scipy BDF/Radau는 수치적으로 동일한 결과(583 FBA calls, biomass=10.01)를 Python 네이티브로 제공
- 공간 시뮬레이션이 현재 연구 범위 외 (well-mixed 가정)
- 향후 공간 시뮬레이션 필요 시 COMETS로 전환 가능

### 3. dFBA 수치 안정성 — Adaptive time-stepping, stiff ODE 해법

**검증 결과:**

| Solver | Time(s) | FBA_calls | Final_biomass | 판정 |
|--------|---------|-----------|---------------|------|
| **BDF** | 0.78 | 583 | 10.01 | 정확 (기본 권장) |
| Radau | 1.69 | 1,504 | 10.01 | 정확 (느림) |
| Euler | 0.16 | 50 | **57.34 (5.7x 과대)** | **사용 금지** |

**Stiff ODE 해법 적용 사례:**
- dFBA는 biomass dynamics가 exponential phase에서 급격히 변화 → stiff system
- scipy.integrate.solve_ivp의 BDF(Backward Differentiation Formula)는 stiff system에 표준적 선택
- Radau(Implicit Runge-Kutta)도 stiff system에 적합하나 BDF 대비 2-3배 느림
- Adaptive time-stepping: solve_ivp가 기본 지원 (rtol=1e-6, atol=1e-8 기본값, 조정 가능)

**COMETS의 수치 방법 (웹 페칭 확인):**
- "discrete approximation of diffusion" — 공간 확산의 이산 근사
- 각 time step에서 dFBA로 개별 종의 대사 활동을 계산
- 자체 adaptive time-stepping 구현 추정 (상세는 Harcombe 2014 참조)

**설계 결정**: BDF 기본, Radau 보조 검증
- 근거: BDF가 Radau 대비 2배 빠르면서 동일 정확도
- Euler는 절대 사용 금지 (5.7x 과대추정 관측)

### 4. NSGA-II + 미생물 군집 — 접종 비율 최적화 적용 사례

**pymoo 라이브러리 (2026-04-26 웹 페칭 확인):**
- NSGA-II 구현: pymoo (Apache-2.0 라이선스)
- 논문: "J. Blank and K. Deb, pymoo: Multi-Objective Optimization in Python, IEEE Access, 2020"
- 기능:
  - non-dominated sorting으로 Pareto fronts 구성
  - crowding distance로 다양성 유지
  - genetic operators (selection, crossover, mutation)로 해 공간 탐색
  - 결과는 Pareto-optimal solution set (단일 해가 아닌 trade-off 옵션 집합)
- 기본 설정: pop_size=100, 200 generations → Pareto front 30해 도출 확인
- 문서: https://pymoo.org/

**미생물 군집 접종 비율 최적화 적용:**
- NSGA-II의 다목적 최적화는 서로 경쟁하는 목적(biomass 생산 vs 대사산물 수율 등)의 Pareto front를 도출
- 접종 비율(inoculation ratio)은 연속 변수로 정의, 각 종의 초기 비율이 최적화 변수
- pymoo의 Problem 클래스를 상속하여 dFBA 시뮬레이터를 objective function으로 래핑

**계산비용 리스크:**
- NSGA-II 1세대 = pop_size(100) × objective evaluations
- 각 evaluation = 1회 dFBA 시뮬레이션 (BDF: ~0.78s)
- 200세대 × 100개체 = 20,000 evaluations × 0.78s = ~4.3h
- 완화: pop_size=30, generations=50으로 축소 → 1,500 × 0.78s = ~20min (초기 탐색)

**관련 연구 사례 (문헌 기반):**
- NSGA-II는 대사 공학에서 유전자 knockout 최적화에 널리 적용
- 접종 비율 최적화에 직접 적용된 사례는 제한적이나, pymoo의 Problem 클래스 래핑으로 구현 가능
- Surrogate-assisted NSGA-II: Phase 3의 surrogate model로 dFBA evaluation 일부 대체 → 계산비용 70-90% 절감 기대

---

### Phase 4 종합 평가

**기술적 실현성 등급: MEDIUM**

**핵심 리스크 3가지:**
1. **NSGA-II + dFBA 계산비용**: 1세대당 100회 dFBA × 0.78s = 78s. 200세대 = 4.3h. 대형 모델(iJO1366)에서는 FBA당 2-5배 느려져 8-20h 예상
2. **dFBA 초기 조건 민감도**: 접종 비율, 초기 농도, pH 등 초기 조건이 시뮬레이션 결과에 큰 영향. 잘못된 초기 조건은 infeasible solution 유발
3. **NSGA-II 수렴 불확실성**: Pareto front의 질이 목적 함수의 복잡도에 따라 달라짐. dFBA objective는 비선형+불연속(조건부 활성/비활성 반응)이어서 수렴 보장 어려움

**대안/완화 전략:**
1. 계산비용 → (a) 초기 탐색은 소규모 pop(30)/gen(50), (b) surrogate model(Module B)로 dFBA 일부 대체, (c) COBRApy multiprocessing으로 병렬 evaluation
2. 초기 조건 → (a) COBRApy 기본 경계 조건을 베이스라인으로 사용, (b) Latin Hypercube Sampling으로 초기 조건 공간 탐색
3. 수렴 → (a) pymoo의 convergence metric(Hypervolume) 모니터링, (b) NSGA-III(참조점 기반) 대안 검토, (c) surrogate-assisted NSGA-II로 evaluation 수 축소

**필요 컴퓨팅 자원:**
- CPU: 8-core 이상 (NSGA-II 병렬 evaluation)
- RAM: 16GB 이상 (dFBA + COBRApy + pymoo 동시 실행)
- GPU: 불필요 (dFBA는 CPU 연산)
- 예상 총시간: 소규모 초기 탐색 ~20min, 전체 실행 ~4-5h (textbook), ~10-20h (iJO1366)

**오픈소스 도구/프레임워크:**
| 도구 | 버전 | URL | 용도 | 상태 |
|------|------|-----|------|------|
| COMETS | v2.12.4 | https://github.com/segrelab/COMETS | 공간 dFBA 시뮬레이션 (Java) | 활성 (2025-06 릴리즈, 380 commits) |
| cometspy | v0.6.1 | https://github.com/segrelab/cometspy | COMETS Python 인터페이스 | 활성 (Virus/Chemotaxis 모델 추가) |
| pymoo | v0.6.1.6 | https://github.com/anyoptimization/pymoo | NSGA-II 다목적 최적화 | 활성 (IEEE Access 2020) |
| scipy BDF/Radau | 1.17.1 | https://scipy.org/ | Stiff ODE solver (dFBA) | 검증 완료 |
| COBRApy | v0.31.1 | https://github.com/opencobra/cobrapy | FBA 실행, 모델 관리 | 활성 (2026-03 릴리즈) |

---

## Phase 3-4 통합 평가

| 항목 | Phase 3 | Phase 4 |
|------|---------|---------|
| 실현성 등급 | MEDIUM | MEDIUM |
| 최대 리스크 | GNN R2 < 0 (데이터+pretraining) | NSGA-II+dFBA 계산비용 |
| 선행 조건 | COBRApy 병렬 FBA 구축 | Phase 3 surrogate 안정화 (R2 > 0.3) |
| 예상 소요 시간 | 1-2일 (초기), 1주 (AL 루프) | 1일 (소규모), 3-5일 (전체) |
| 성공 기준 | Surrogate R2 > 0.5 | Pareto front 30해 + TOPSIS tau > 0.7 |

**의존 관계**: Phase 4의 NSGA-II는 Phase 3의 surrogate model이 안정화된 후에야 surrogate-assisted optimization 가능. Phase 3 실패 시 Phase 4는 brute-force dFBA + NSGA-II로 진행 가능하나 계산비용 급증.

---

## v2 업데이트 내역 (run_01 대비 변경점)

| 항목 | run_01 | run_02 (본 문서) | 변경 근거 |
|------|--------|----------------|----------|
| PyG 이종 GNN 레이어 | 3개 나열 | 4개 나열 + 동종 20+ 레이어 명시 | 웹 페칭으로 PyG 저장소 확인 |
| COBRApy 버전 | v0.31.1 | v0.31.1 (릴리즈일 2026-03-26 명시) | 웹 페칭으로 cobrapy 저장소 확인 |
| COMETS 상세 | 기본 정보 | v2.12.4 (2025-06-18), 코드 구성, cometspy v0.6.1 상세 | 웹 페칭으로 COMETS/cometspy 저장소 확인 |
| FLYCOP | 404 | 404 재확인 + 논문 DOI 페칭도 실패 | 웹 페칭 재시도 3회 모두 실패 |
| dfba-python | FAIL | FAIL 재확인 | 웹 페칭 (PyPI + GitHub) 모두 404 |
| cFBA | 확인 불가 | 확인 불가 재확인 | PyPI JavaScript challenge, GitHub 404 |
| pymoo | 기본 정보 | NSGA-II 동작 상세 (non-dominated sorting, crowding distance, genetic operators) | 웹 페칭으로 pymoo 저장소 확인 |
| Module A | 미구현 | 구현 완료 언급 | outputs/analysis/run_01/module_a_fba_generator.py 작성됨 |
| 도구 테이블 | URL만 | 버전+URL+상태 열 추가 | 웹 페칭으로 실제 버전/상태 확인 |

---

## 웹 검색/페칭 결과 로그

**성공한 페칭:**
1. https://github.com/opencobra/cobrapy → COBRApy v0.31.1, 2026-03-26 릴리즈, 병렬 미지원 확인
2. https://github.com/segrelab/COMETS → COMETS v2.12.4, 2025-06-18, Java 기반, dFBA+diffusion
3. https://github.com/segrelab/cometspy → cometspy v0.6.1, Python 인터페이스, pip3 설치
4. https://github.com/pyg-team/pytorch_geometric → PyG 이종 GNN 4+ 레이어, 동종 20+ 레이어, multi-GPU 지원
5. https://github.com/anyoptimization/pymoo → NSGA-II 상세 동작, IEEE Access 2020 논문
6. https://segrelab.github.io/comets-manual/ → COMETS manual, discrete diffusion, dFBA 시뮬레이션 예시

**실패한 페칭:**
1. https://github.com/baliga-lab/FLYCOP → 404 (저장소 삭제)
2. https://github.com/baliga-lab/Flycop → 404
3. https://github.com/biosustain/dfba → 404 (저장소 삭제)
4. https://pypi.org/project/dfba/ → JavaScript challenge (내용 로드 불가)
5. https://pypi.org/project/cfba/ → JavaScript challenge (내용 로드 불가)
6. https://github.com/araitzelmo/muBialSim → 404
7. https://link.springer.com/article/10.1186/s12918-018-0639-6 → 303 리다이렉트 루프
8. https://pmc.ncbi.nlm.nih.gov/articles/PMC6265957/ → 잘못된 논문 반환 (화학 논문)
9. https://pmc.ncbi.nlm.nih.gov/articles/PMC6201919 → 잘못된 논문 반환 (condensin inhibitors)
10. arXiv 무작위 ID 2개 → 관련 없는 논문 (수학/물리)
11. WebSearch 4회 → "Web search temporarily unavailable" (전체 실패)
12. BrainSearch → 도구 오류 (tool.mapToolResultToToolResultBlockParam is not a function)

**결론**: WebSearch가 완전히 불가능하여 학술 논문 직접 검색은 수행하지 못함. 그러나 GitHub 저장소 페칭으로 핵심 도구들의 버전, 기능, 활성도는 확인 완료. 문헌 검색은 향후 WebSearch 복구 후 보강 필요.
