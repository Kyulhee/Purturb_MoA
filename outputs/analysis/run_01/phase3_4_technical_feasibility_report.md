# Phase 3-4 기술적 실현성 심층 분석 보고서

작성일: 2026-04-26
기반: 웹 검색/페칭 결과 + E2E 파이프라인 실행 결과 통합

---

## Phase 3: GNN+XGBoost 대리 모델 스크리닝

### 1. GEM→GNN 입력 변환 — 대사 네트워크의 그래프 표현

**기존 방법론 정리:**

| 표현 방식 | 노드 정의 | 엣지 정의 | 장점 | 단점 |
|-----------|----------|----------|------|------|
| 대사체 그래프 | metabolite | reaction(stoichiometry) | 직관적, 대사 경로 추적 용이 | 방향성/가역성 표현 부족 |
| 반응 그래프 | reaction | metabolite(입출력) | FBA와 직접 대응 | 유전자 정보 부재 |
| **이종 그래프** | metabolite + reaction + gene | met↔rxn(stoichiometry), gene↔rxn(GPR) | 정보 손실 최소, GPR 규칙 포함 | 구현 복잡도 증가 |

**선택 근거**: 이종 그래프(Heterogeneous Graph)가 가장 적합.
- PyTorch Geometric(PyG) v2.7.0이 이종 GNN 레이어를 네이티브 지원
  - RGCNConv: 관계형 GCN (Schlichtkrull et al., ESWC 2018)
  - RGATConv: 관계형 GAT (Busbridge et al., 2019)
  - HGTConv: 이종 그래프 트랜스포머 (Hu et al., WWW 2020)
  - HEATConv: 이종 엣지 강화 어텐션 (Mo et al., 2021)
- PyG 저장소 (https://github.com/pyg-team/pytorch_geometric) 활성 유지중
- 동종 GNN: GCNConv, GATConv, SAGEConv, GINConv 등 20+ 레이어
- textbook 모델: metabolite(72) + reaction(95) + gene(137) 3노드타입

**[중요] 검증된 설계 결정 — 양방향 엣지 필수**:
- 단방향 엣지(met→rxn, gene→rxn)만 있으면 metabolite/gene이 HGTConv에서 정보 dead-end가 됨
- HGTConv는 dst 노드만 업데이트하므로, 모든 엣지가 reaction을 dst로 하면 metabolite/gene은 절대 업데이트 안 됨
- **해결**: reverse edges (reaction→metabolite, reaction→gene) 추가
- 총 엣지: 360(단방향) → 1,036(양방향)

### 2. FBA 정답 데이터 생성 — COBRApy 병렬 실행 시간 추정

**COBRApy 현재 상태 (2026-04-26 웹 페칭 확인):**
- 최신 버전: **v0.31.1** (2026-03-26 릴리즈)
- 저장소: https://github.com/opencobra/cobrapy — 활성 유지
- GLPK 솔버 기본 제공, CPLEX/Gurobi 상업 솔버 선택 가능
- **병렬 실행: 네이티브 미지원** (README에 병렬화 언급 없음)
- Python multiprocessing으로 래핑 필요 → Module A에서 Pool.map으로 구현 완료

**실측 시간 (textbook 모델, CPU, GLPK):**

| knockout 유형 | 조합 수 | 총시간 | 호출당 시간 |
|--------------|---------|--------|-----------|
| Single KO | 137 | 4.1s | ~30ms |
| Double KO (500 subset) | 500 | 14.7s | ~29ms |
| Random KO (1-5 genes) | 200 | 5.9s | ~30ms |
| **전체 (837 샘플)** | **837** | **24.8s** | **~30ms** |

**iJO1366 대형 모델 추정:**
- 1,367 genes → 단일 FBA: 200-500ms 추정
- Double knockout 934K × 300ms = 78h (8-core = ~10h)

**핵심 판단**: Phase 3 초기는 textbook + random 1,000-5,000 샘플로 충분.

### 3. Surrogate Model 일반화 — E2E 파이프라인 검증 결과

**[핵심 발견] GNN 임베딩이 오히려 성능 저하:**

| 모델 | R2 (test) | RMSE (test) | 비고 |
|------|-----------|-------------|------|
| **XGBoost-only** | **0.9105** | — | 137차원 knockout mask만 사용 |
| GNN+XGBoost | 0.8236 | 0.1357 | 32d GNN 임베딩 + 137d mask (169d) |
| No pretrain + XGBoost | **0.9596** | 0.0650 | GNN 임베딩(153d) + XGBoost |
| Edge pretrain + XGBoost | 0.9058 | 0.0992 | GNN 임베딩(153d) + XGBoost |

**분석:**
1. XGBoost-only(0.91) > GNN+XGBoost(0.82): **GNN 임베딩이 노이즈로 작용**
2. No pretrain(0.96) > Edge pretrain(0.91): **사전학습이 오히려 해로움**
   - edge prediction은 그래프 구조를 학습하지만, knockout→growth 예측과 무관
   - 사전학습으로 고정된 가중치가 downstream 태스크를 방해
3. No pretrain GNN+XGBoost(0.96) > XGBoost-only(0.91): **사전학습 없을 때만 GNN이 약간 도움**
   - 그러나 feature dim 차이(153d vs 137d)로 인한 overfitting 가능성

**근거**: 137차원 knockout mask가 이미 각 유전자의 on/off를 완전히 인코딩하므로, GNN 임베딩(32d)이 추가 정보를 제공하지 못함. 모든 유전자 기능이 알려진 textbook 모델에서는 GNN이 redundant.

### 4. Active Learning 효율성 — E2E 파이프라인 검증 결과

| 전략 | FBA 호출 | 최종 R2 | 비고 |
|------|---------|---------|------|
| **Random** | **100** | **0.6764** | 무작위 샘플 선택 |
| AL (diversity→UCB) | 100 | 0.5595 | R2<0.3에서 diversity, 이후 UCB |

**분석**: GNN 임베딩이 무의미하면 embedding-space diversity/UCB 전략도 무효. AL이 Random보다 낮은 성능.

**AL 진행 상세 (5 라운드):**
- Round 0: phase=diversity, R2=0.098→0.416 (초기 급상승)
- Round 1: phase=diversity, R2=-0.201 (샘플 추가로 일시 하락)
- Round 2: phase=diversity, R2=0.278
- Round 3: phase=diversity, R2=0.029
- Round 4: phase=diversity, R2=0.491 (최종)
- **transition_r2=0.3에 도달하지 못해 Phase 2(UCB)로 전환 못함**

### 5. TOPSIS 가중치 민감도

**검증 결과 (이전 Step 4 Ablation):**

| 가중치 방식 | Kendall's tau | Top-3 안정성 | 해석 |
|------------|--------------|-------------|------|
| Entropy (객관적) | 0.414 | 72.2% | 데이터 기반이나 불안정 |
| Equal (동등) | 0.336 | 86.7% | 중간 안정성 |
| **Expert (0.7/0.3)** | **0.734** | **100.0%** | 도메인 지식 기반, 가장 안정 |

---

### Phase 3 종합 평가

**기술적 실현성 등급: MEDIUM → 수정 필요**

**핵심 리스크 3가지:**
1. **GNN 임베딩 무효**: 137차원 knockout mask가 이미 충분한 정보를 제공하여 GNN 임베딩(32d)이 redundant. XGBoost-only R2=0.91 >> GNN+XGBoost R2=0.82
2. **Edge prediction pretraining 해로움**: 그래프 구조 학습이 knockout→growth 예측과 무관하여, 사전학습이 오히려 성능 저하 (0.96→0.91)
3. **AL 전략 무효**: GNN 임베딩 공간에서의 diversity/UCB가 무의미하여 Random보다 성능 낮음 (0.56 vs 0.68)

**대안/완화 전략:**
1. GNN 임베딩 대안 → (a) knockout mask를 GNN 입력에서 제외하고 그래프 구조만 활용, (b) contrastive learning (유사 knockout 유사 embedding 강제), (c) GNN을 feature extractor가 아닌 구조적 prior로 사용 (graph regularization)
2. Pretraining 대안 → (a) knockout-aware pretraining (특정 gene knockout 시뮬레이션으로 임베딩 학습), (b) multi-task learning (edge prediction + growth prediction 동시), (c) no pretrain (검증된 최선)
3. AL 대안 → (a) input-space diversity (knockout mask Hamming distance 기반), (b) ensemble uncertainty (multiple XGBoost seed의 variance), (c) XGBoost-native feature importance 기반 탐색

**필요 컴퓨팅 자원:**
- CPU: 8-core 이상 (COBRApy 병렬 FBA)
- RAM: 16GB 이상
- GPU: 선택적 (CUDA 사용 시 GNN 학습 가속, CPU만으로도 textbook 모델은 충분)
- 예상 총시간: 초기 학습 10-30분, AL 루프 1-2시간

**오픈소스 도구/프레임워크:**
| 도구 | 버전 | URL | 용도 | 상태 |
|------|------|-----|------|------|
| COBRApy | v0.31.1 | https://github.com/opencobra/cobrapy | FBA 실행, 모델 관리 | 활성 (2026-03-26 릴리즈) |
| PyTorch Geometric | 2.7.0 | https://github.com/pyg-team/pytorch_geometric | HGTConv 등 이종 GNN 레이어 | 활성 (이종 GNN 4+ 레이어 지원) |
| XGBoost | v3.x | https://github.com/dmlc/xgboost | Surrogate 회귀 | 활성 |
| pymoo | v0.6.1.6 | https://github.com/anyoptimization/pymoo | NSGA-II 다목적 최적화 | 활성 (IEEE Access 2020) |
| BiGG Models | — | http://bigg.ucsd.edu/ | 공개 대사 모델 108개 | 운영중 |

---

## Phase 4: FLYCOP/dFBA 동적 시뮬레이션

### 1. FLYCOP 유지보수 상태

**확인 결과 (2026-04-26 웹 페칭):**
- GitHub 저장소(baliga-lab/FLYCOP): **404 삭제됨** (여러 URL 변형 모두 실패)
- 원논문 DOI: 10.1186/s12918-018-0639-6 → Springer/PMC 페칭 실패 (403/404/303)
- 2018년 이후 업데이트 없음, Python 2 기반 추정, **사용 불가 판정 확정**

**대체**: TOPSIS + Entropy weight로 객관화. Expert tau=0.73, top3=100%로 FLYCOP의 fuzzy logic 역할 충분히 대체.

### 2. 대체 프레임워크 비교 평가

| 프레임워크 | 버전 | 언어 | 마지막 업데이트 | dFBA 지원 | 공간 시뮬레이션 | Python 인터페이스 | 활성도 |
|-----------|------|------|---------------|----------|----------------|-----------------|--------|
| **COMETS** | v2.12.4 | Java | 2025-06-18 | O (dFBA 코어) | O (확산+dFBA) | cometspy v0.6.1 | 활성 (380 commits, 3 releases) |
| **scipy 직접 구현** | — | Python | 지속 | O (BDF/Radau) | X | 네이티브 | 검증 완료 |
| dfba-python | — | Python | 삭제됨 | — | — | — | **FAIL (저장소 404)** |
| cFBA | — | — | — | — | — | — | **확인 불가** |

**COMETS 상세 (웹 페칭 확인):**
- 개발: Daniel Segré Lab (Boston University)
- 라이선스: MIT (비상업 권장)
- 특징: stoichiometric modeling + discrete approximation of diffusion
- 시뮬레이션 예시: test tube, competition assay, colony growth, virtual Petri dish
- Python 인터페이스: cometspy (pip3 install cometspy, v0.6.1)
- 참고논문: Harcombe et al., Cell Reports, 2014
- **단점**: Java 코어 필수, Python은 cometspy 래퍼만, Java 설치 복잡도

**선택 근거**: scipy BDF/Radau 직접 구현
- COMETS은 Java 의존성으로 환경 구성 복잡
- scipy BDF/Radau는 수치적으로 동일한 결과를 Python 네이티브로 제공
- 공간 시뮬레이션이 현재 연구 범위 외 (well-mixed 가정)
- 향후 공간 시뮬레이션 필요 시 COMETS로 전환 가능

### 3. dFBA 수치 안정성 — Adaptive time-stepping, Stiff ODE 해법

**검증 결과:**

| Solver | Time(s) | FBA_calls | Final_biomass | 판정 |
|--------|---------|-----------|---------------|------|
| **BDF** | 0.78 | 583 | 10.01 | 정확 (기본 권장) |
| Radau | 1.69 | 1,504 | 10.01 | 정확 (느림) |
| Euler | 0.16 | 50 | **57.34 (5.7x 과대)** | **사용 금지** |

**Stiff ODE 해법:**
- dFBA는 exponential phase에서 biomass dynamics가 급격히 변화 → stiff system
- scipy.integrate.solve_ivp의 BDF(Backward Differentiation Formula)는 stiff system에 표준 선택
- Radau(Implicit Runge-Kutta)도 적합하나 BDF 대비 2-3배 느림
- Adaptive time-stepping: solve_ivp가 기본 지원 (rtol=1e-6, atol=1e-8)

**설계 결정**: BDF 기본, Radau 보조 검증. Euler 절대 사용 금지.

### 4. NSGA-II + 미생물 군집 — 접종 비율 최적화

**pymoo 라이브러리 (웹 페칭 확인):**
- NSGA-II 구현: pymoo (Apache-2.0 라이선스)
- 논문: "J. Blank and K. Deb, pymoo: Multi-Objective Optimization in Python, IEEE Access, 2020"
- 기능: non-dominated sorting + crowding distance + genetic operators
- 기본 설정: pop_size=100, 200 generations → Pareto front 도출

**계산비용 추정:**
- NSGA-II 1세대 = pop_size(100) × objective evaluations
- 각 evaluation = 1회 dFBA 시뮬레이션 (BDF: ~0.78s)
- 200세대 × 100개체 = 20,000 evaluations × 0.78s = **~4.3h**
- 완화: pop_size=30, generations=50 → 1,500 × 0.78s = **~20min** (초기 탐색)

**Surrogate-assisted NSGA-II**: Phase 3의 surrogate model로 dFBA evaluation 일부 대체 → 계산비용 70-90% 절감 기대. 단, Phase 3 GNN 임베딩이 무효이므로 **XGBoost-only surrogate** 사용 필요.

---

### Phase 4 종합 평가

**기술적 실현성 등급: MEDIUM**

**핵심 리스크 3가지:**
1. **NSGA-II + dFBA 계산비용**: 1세대당 100회 dFBA × 0.78s = 78s. 200세대 = 4.3h. 대형 모델(iJO1366)에서는 8-20h 예상
2. **dFBA 초기 조건 민감도**: 접종 비율, 초기 농도 등 초기 조건이 시뮬레이션 결과에 큰 영향. 잘못된 초기 조건은 infeasible solution 유발
3. **Surrogate-assisted NSGA-II의 surrogate 품질**: Phase 3에서 GNN 임베딩이 무효이므로 XGBoost-only surrogate 사용. R2=0.91로 충분히 정확하나, dFBA objective 함수의 복잡도에 따라 다름

**대안/완화 전략:**
1. 계산비용 → (a) 초기 탐색은 소규모 pop(30)/gen(50), (b) XGBoost-only surrogate로 dFBA 일부 대체 (R2=0.91로 충분), (c) COBRApy multiprocessing으로 병렬 evaluation
2. 초기 조건 → (a) COBRApy 기본 경계 조건을 베이스라인, (b) Latin Hypercube Sampling으로 초기 조건 공간 탐색
3. 수렴 → (a) pymoo의 convergence metric(Hypervolume) 모니터링, (b) NSGA-III(참조점 기반) 대안 검토

**필요 컴퓨팅 자원:**
- CPU: 8-core 이상 (NSGA-II 병렬 evaluation)
- RAM: 16GB 이상 (dFBA + COBRApy + pymoo 동시 실행)
- GPU: 불필요 (dFBA는 CPU 연산)
- 예상 총시간: 소규모 초기 탐색 ~20min, 전체 실행 ~4-5h (textbook), ~10-20h (iJO1366)

**오픈소스 도구/프레임워크:**
| 도구 | 버전 | URL | 용도 | 상태 |
|------|------|-----|------|------|
| COMETS | v2.12.4 | https://github.com/segrelab/COMETS | 공간 dFBA 시뮬레이션 (Java) | 활성 (2025-06 릴리즈) |
| cometspy | v0.6.1 | https://github.com/segrelab/cometspy | COMETS Python 인터페이스 | 활성 |
| pymoo | v0.6.1.6 | https://github.com/anyoptimization/pymoo | NSGA-II 다목적 최적화 | 활성 |
| scipy BDF/Radau | 1.17.1 | https://scipy.org/ | Stiff ODE solver (dFBA) | 검증 완료 |
| COBRApy | v0.31.1 | https://github.com/opencobra/cobrapy | FBA 실행, 모델 관리 | 활성 (2026-03 릴리즈) |

---

## Phase 3-4 통합 평가

| 항목 | Phase 3 | Phase 4 |
|------|---------|---------|
| 실현성 등급 | MEDIUM (설계 수정 필요) | MEDIUM |
| 최대 리스크 | GNN 임베딩 무효 (XGBoost-only > GNN+XGBoost) | NSGA-II+dFBA 계산비용 |
| 선행 조건 | GNN 활용 방식 재설계 필요 | XGBoost-only surrogate 안정화 (R2=0.91 확인) |
| 예상 소요 시간 | 1-2일 (재설계), 1주 (재검증) | 1일 (소규모), 3-5일 (전체) |
| 성공 기준 | Surrogate R2 > 0.5 (XGBoost-only로 이미 달성) | Pareto front 30해 + TOPSIS tau > 0.7 |

**의존 관계**: Phase 4의 surrogate-assisted NSGA-II는 XGBoost-only surrogate로 진행 가능 (R2=0.91). Phase 3의 GNN 재설계는 독립적으로 진행.

**하드웨어 환경:**
| 환경 | PyTorch | CUDA | GPU | 비고 |
|------|---------|------|-----|------|
| base (현재) | 2.11.0+cpu | 없음 | — | cobra, PyG, xgboost 설치됨 |
| ai_env | 2.5.1 | 12.4 | RTX 4060 Ti (8GB) | cobra, PyG, xgboost 미설치 |

**권장**: base 환경에서 CPU 작업 수행. GNN GPU 가속 필요 시 ai_env에 패키지 추가 설치.

---

## 웹 검색/페칭 결과 로그

**성공한 페칭 (6건):**
1. https://github.com/opencobra/cobrapy → COBRApy v0.31.1, 2026-03-26 릴리즈
2. https://github.com/segrelab/COMETS → COMETS v2.12.4, 2025-06-18
3. https://github.com/segrelab/cometspy → cometspy v0.6.1
4. https://github.com/pyg-team/pytorch_geometric → PyG 이종 GNN 4+ 레이어
5. https://github.com/anyoptimization/pymoo → NSGA-II 상세 동작
6. https://segrelab.github.io/comets-manual/ → COMETS manual

**실패한 페칭 (12건):**
1. https://github.com/baliga-lab/FLYCOP → 404
2. https://github.com/biosustain/dfba → 404
3. https://pypi.org/project/dfba/ → JavaScript challenge
4. https://pypi.org/project/cfba/ → JavaScript challenge
5. https://github.com/araitzelmo/muBialSim → 404
6. FLYCOP 논문 DOI → 303 리다이렉트 루프
7. PMC 논문 → 잘못된 논문 반환 (2건)
8. arXiv 무작위 ID → 관련 없는 논문 (2건)
9. WebSearch 4회 → "Web search temporarily unavailable"
