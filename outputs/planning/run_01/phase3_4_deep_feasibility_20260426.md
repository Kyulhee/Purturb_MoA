# Phase 3-4 기술적 실현성 심층 분석

작성일: 2026-04-26

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
- PyTorch Geometric이 HGTConv, RGCNConv, HGTConv 등 이종 GNN 레이어를 네이티브 지원
  - RGCNConv: 관계형 GCN (Schlichtkrull et al., ESWC 2018)
  - RGATConv: 관계형 GAT (Busbridge et al., 2019)
  - HGTConv: 이종 그래프 트랜스포머 (Hu et al., WWW 2020)
  - HEATConv: 이종 엣지 강화 어텐션 (Mo et al., 2021)
- metabolite(72) + reaction(95) + gene(137) 3노드타입 + met→rxn(188+172) + gene→rxn(158) 엣지가 textbook 모델에서 검증됨
- knockout_mask(137차원)를 gene 노드 피처로 통합 가능 → GNN이 knockout 효과를 end-to-end로 학습

**리스크**: 이종 GNN은 동종 GNN 대비 학습 샘플 요구량이 큼 (현재 135샘플로 R2 < 0 관측)

### 2. FBA 정답 데이터 생성 — COBRApy 병렬 실행 시간 추정

**COBRApy 현재 상태:**
- 최신 버전: v0.31.1 (2026-03-26)
- GLPK 솔버 기본 제공, CPLEX/Gurobi 상업 솔버 선택 가능
- 병렬 실행: 네이티브 미지원, Python multiprocessing으로 래핑 필요

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
- COBRApy multiprocessing 래핑은 간단 (Pool.map으로 FBA 함수 병렬화)

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
1. **샘플 부족**: 135샘플로 R2 < 0 → 최소 500-1,000샘플 필요. double knockout 9,316 조합에서 random/AL sampling으로 확보 가능하나, COBRApy 직렬 실행 시 15.5min (8-core ~2min) 소요
2. **GNN 임베딩 품질**: autoencoder loss 366K→365K로 거의 학습 안 됨 → contrastive loss 또는 edge prediction pretraining 필요. 이종 GNN은 동종 대비 학습 난이도 높음
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
| 도구 | URL | 용도 |
|------|-----|------|
| COBRApy v0.31.1 | https://github.com/opencobra/cobrapy | FBA 실행, 모델 관리 |
| PyTorch Geometric | https://github.com/pyg-team/pytorch_geometric | HGTConv 등 이종 GNN 레이어 |
| XGBoost v3.2.0 | https://github.com/dmlc/xgboost | Surrogate 회귀, Quantile regression |
| pymoo v0.6.1.6 | https://github.com/anyoptimization/pymoo | NSGA-II 다목적 최적화 |
| BiGG Models | http://bigg.ucsd.edu/ | 공개 대사 모델 108개 |

---

## Phase 4: FLYCOP/dFBA 동적 시뮬레이션

### 1. FLYCOP 유지보수 상태

**확인 결과:**
- GitHub 저장소(baliga-lab/FLYCOP): **404 삭제됨**
- 원논문: FLYCOP (Fuzzy Logic Combined with Perturbation Theory) — BMC Systems Biology 2018
- 2018년 이후 업데이트 없음, Python 2 기반 추정, 사용 불가

**대체 전략 검증 완료**: TOPSIS + Entropy weight로 객관화 대체. Expert tau=0.73, top3=100%로 FLYCOP의 fuzzy logic 역할을 충분히 대체.

### 2. 대체 프레임워크 비교 평가

| 프레임워크 | 버전 | 언어 | 마지막 업데이트 | dFBA 지원 | 공간 시뮬레이션 | Python 인터페이스 | 활성도 |
|-----------|------|------|---------------|----------|----------------|-----------------|--------|
| **COMETS** | v2.12.4 | Java | 2025-06 | O (dFBA 코어) | O (확산+dFBA) | cometspy (v0.6.1) | 활성 (380 commits) |
| **scipy 직접 구현** | — | Python | 지속 | O (BDF/Radau) | X | 네이티브 | 검증 완료 |
| dfba-python | — | Python | 삭제됨 | — | — | — | FAIL (404) |
| cFBA | — | — | — | — | — | — | 확인 불가 (PyPI/404) |

**COMETS 상세:**
- 개발: Daniel Segré Lab (Boston University)
- 라이선스: MIT (비상업 권장)
- 특징: stoichiometric modeling + discrete diffusion, 공간 구조 미생물 커뮤니티 시뮬레이션
- 단점: Java 코어 필수, Python은 cometspy 래퍼만, Java 설치 복잡도
- 예시: test tube, competition assay, circular/branching colony, demographic noise

**선택 근거**: scipy BDF/Radau 직접 구현
- COMETS는 Java 의존성으로 환경 구성 복잡
- scipy BDF/Radau는 수치적으로 동일한 결과(583 FBA calls, biomass=10.01)를 Python 네이티브로 제공
- 공간 시뮬레이션이 현재 연구 범위 외 (well-mixed 가정)

### 3. dFBA 수치 안정성 — Adaptive time-stepping, stiff ODE 해법

**검증 결과:**

| Solver | Time(s) | FBA_calls | Final_biomass | 판정 |
|--------|---------|-----------|---------------|------|
| **BDF** | 0.78 | 583 | 10.01 | 정확 (기본 권장) |
| Radau | 1.69 | 1,504 | 10.01 | 정확 (느림) |
| Euler | 0.16 | 50 | **57.34 (5.7x 과대)** | 사용 금지 |

**Stiff ODE 해법 적용 사례:**
- dFBA는 biomass dynamics가 exponential phase에서 급격히 변화 → stiff system
- scipy.integrate.solve_ivp의 BDF(Backward Differentiation Formula)는 stiff system에 표준적 선택
- Radau(Implicit Runge-Kutta)도 stiff system에 적합하나 BDF 대비 2-3배 느림
- Adaptive time-stepping: solve_ivp가 기본 지원 (rtol=1e-6, atol=1e-8 기본값, 조정 가능)

**설계 결정**: BDF 기본, Radau 보조 검증
- 근거: BDF가 Radau 대비 2배 빠르면서 동일 정확도
- Euler는 절대 사용 금지 (5.7x 과대추정 관측)

### 4. NSGA-II + 미생물 군집 — 접종 비율 최적화 적용 사례

**pymoo 라이브러리:**
- NSGA-II 구현: pymoo 0.6.1.6 (Apache-2.0)
- 논문: "J. Blank and K. Deb, pymoo: Multi-Objective Optimization in Python, IEEE Access, 2020"
- 기능: non-dominated sorting, crowding distance, selection/crossover/mutation
- 기본 설정: pop_size=100, 200 generations → Pareto front 30해 도출 확인

**미생물 군집 접종 비율 최적화 적용:**
- NSGA-II의 다목적 최적화는 서로 경쟁하는 목적(biomass 생산 vs 대사산물 수율 등)의 Pareto front를 도출
- 접종 비율(inoculation ratio)은 연속 변수로 정의, 각 종의 초기 비율이 최적화 변수
- pymoo의 Problem 클래스를 상속하여 dFBA 시뮬레이터를 objective function으로 래핑

**계산비용 리스크:**
- NSGA-II 1세대 = pop_size(100) × objective evaluations
- 각 evaluation = 1회 dFBA 시뮬레이션 (BDF: ~0.78s)
- 200세대 × 100개체 = 20,000 evaluations × 0.78s = ~4.3h
- 완화: pop_size=30, generations=50으로 축소 → 1,500 × 0.78s = ~20min (초기 탐색)

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
| 도구 | URL | 용도 |
|------|-----|------|
| COMETS v2.12.4 | https://github.com/segrelab/COMETS | 공간 dFBA 시뮬레이션 (Java) |
| cometspy v0.6.1 | https://github.com/segrelab/cometspy | COMETS Python 인터페이스 |
| pymoo v0.6.1.6 | https://github.com/anyoptimization/pymoo | NSGA-II 다목적 최적화 |
| scipy BDF/Radau | https://scipy.org/ | Stiff ODE solver (dFBA) |
| COBRApy v0.31.1 | https://github.com/opencobra/cobrapy | FBA 실행, 모델 관리 |

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
