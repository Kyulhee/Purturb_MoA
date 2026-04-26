# Phase 3-4 통합 아키텍처 상세 설계

> 작성일: 2026-04-26
> 기반: outputs/planning/run_02/phase3_4_feasibility_analysis.md

---

## 1. 전체 아키텍처 개요

```
┌─────────────────────────────────────────────────────────────────┐
│                    OCT LLM XAI Pipeline                         │
│                                                                 │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────────────┐ │
│  │  Phase 3    │    │  Phase 3-4  │    │      Phase 4        │ │
│  │  Surrogate  │◄──►│  Active     │◄──►│  Dynamic Sim        │ │
│  │  Model      │    │  Learning   │    │  (COMETS dFBA)      │ │
│  │  (GNN+XGB)  │    │  Loop       │    │  + NSGA-II          │ │
│  └──────┬──────┘    └──────┬──────┘    └──────────┬──────────┘ │
│         │                  │                      │            │
│  ┌──────▼──────┐    ┌──────▼──────┐    ┌──────────▼──────────┐ │
│  │  FBA Ground │    │  Uncertainty│    │  TOPSIS/Pareto      │ │
│  │  Truth Gen  │    │  Sampling   │    │  Decision Support   │ │
│  │  (COBRApy)  │    │             │    │                     │ │
│  └─────────────┘    └─────────────┘    └─────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
```

---

## 2. 모듈별 상세 설계

### 2.1 모듈 A: FBA Ground Truth Generator

**목적**: 대리 모델 학습용 정답 데이터 대량 생성

**입력**:
- GEM 모델 (E. coli iML1515 또는 대상 미생물 모델)
- 파라미터 공간 정의 (유전자 녹아웃 조합, 환경 조건, 접종 비율)

**출력**:
- (파라미터 조합, FBA 해) 쌍의 데이터셋
- flux distribution, growth rate, byproduct secretion rate

**구현 설계**:

```python
# 의사코드
class FBAGroundTruthGenerator:
    def __init__(self, model_path, solver='glpk', n_workers=8):
        self.model = cobra.io.load_model(model_path)
        self.solver = solver
        self.n_workers = n_workers
    
    def generate_knockout_combinations(self, n_genes, n_combinations, 
                                        max_knockouts=5):
        """유전자 녹아웃 조합 생성"""
        genes = list(self.model.genes)
        combos = []
        for _ in range(n_combinations):
            k = random.randint(1, max_knockouts)
            combo = random.sample(genes, min(k, len(genes)))
            combos.append(combo)
        return combos
    
    def run_fba_single(self, knockout_genes):
        """단일 FBA 실행"""
        with self.model:
            for gene in knockout_genes:
                gene.knock_out()
            solution = self.model.optimize()
            return {
                'knockout_genes': [g.id for g in knockout_genes],
                'growth_rate': solution.objective_value,
                'flux_distribution': solution.fluxes.to_dict(),
                'status': solution.status
            }
    
    def run_parallel(self, combinations, batch_size=1000):
        """병렬 FBA 실행"""
        with multiprocessing.Pool(self.n_workers) as pool:
            results = pool.map(self.run_fba_single, combinations, 
                              chunksize=batch_size // self.n_workers)
        return results
```

**성능 추정** (이전 분석 기반):
| 조합 수 | 8코어 GLPK | 8코어 CPLEX |
|---------|-----------|------------|
| 1,000 | ~5-10분 | ~1-5분 |
| 10,000 | ~30-60분 | ~10-30분 |
| 100,000 | ~5-10시간 | ~2-5시간 |

**데이터 저장 형식**:
- Parquet 형식 (flux vector 저장에 효율적)
- 메타데이터: 모델 ID, 솔버, 타임스탬프

---

### 2.2 모듈 B: GNN+XGBoost Surrogate Model

**목적**: FBA 결과를 빠르게 예측하는 대리 모델

**입력**:
- GEM 모델의 그래프 표현 (HeteroData)
- 유전자 녹아웃/환경 조건 (테이블 데이터)

**출력**:
- 예측 growth rate
- 예측 flux distribution (핵심 flux만)
- 예측 불확실성 (Active Learning용)

**아키텍처**:

```
입력: GEM 그래프 + 노크아웃 마스크
         │
    ┌────▼────┐
    │  GNN    │  ← HGTConv 3층 (이종 그래프)
    │ Encoder │     metabolite/reaction/gene 노드
    └────┬────┘
         │
    graph-level embedding (readout)
         │
    ┌────▼────┐
    │ Concat  │  ← GNN 임베딩 + 테이블 피처(온도, pH, 접종비율 등)
    └────┬────┘
         │
    ┌────▼────┐
    │ XGBoost │  ← 최종 예측 (growth rate, 핵심 flux)
    │ Regressor│     uncertainty: quantile regression
    └─────────┘
```

**GNN 그래프 구성**:

| 노드 타입 | 피처 | 차원 |
|-----------|------|------|
| metabolite | chemical properties (MW, charge, formula encoding) | 32-64 |
| reaction | stoichiometry summary, pathway annotation | 32-64 |
| gene | gene expression level (knockout=0, wt=1) | 1-8 |

| 엣지 타입 | 속성 |
|-----------|------|
| metabolite→reaction | stoichiometric coefficient (-값=substrate) |
| reaction→metabolite | stoichiometric coefficient (+값=product) |
| gene→reaction | GPR rule (AND/OR 인코딩) |

**학습 전략**:
1. Phase 1: GNN 사전학습 (노드 마스킹, 엣지 예측)
2. Phase 2: GNN+XGBoost end-to-end 파인튜닝
3. Phase 3: Active Learning 루프로 반복 개선

**하이퍼파라미터**:
- GNN: HGTConv, 3 layers, hidden_dim=128, heads=4
- XGBoost: n_estimators=500, max_depth=8, learning_rate=0.05
- 학습: AdamW(lr=1e-3), batch_size=32, epochs=100

---

### 2.3 모듈 C: Active Learning Loop (Phase 3-4 연결)

**목적**: 대리 모델과 dFBA 시뮬레이션 간의 효율적인 데이터 수집

**핵심 알고리즘**:

```
초기: 랜덤 샘플링 N₀=1000개 FBA 실행 → 대리 모델 초기 학습

Loop (최대 T회):
  1. 현재 대리 모델로 파라미터 공간 전체 예측 + 불확실성 추정
  2. 불확실성 상위 K=100개 파라미터 조합 선택
     (acquisition function: UCB, EI, 또는 Thompson sampling)
  3. 선택된 K개 조합에 대해 실제 FBA 실행 (ground truth)
  4. 새 데이터로 대리 모델 재학습 (incremental update)
  5. 검증: hold-out set에서 R² > 0.95 달성 시 종료
  6. 미달성 시 NSGA-II 파라미터 공간에서 미탐색 영역 식별 → 2번으로
```

**Acquisition Function 선택**:
- **UCB (Upper Confidence Bound)**: 간단, 탐색-활용 균형
- **EI (Expected Improvement)**: Bayesian Optimization 표준
- **Thompson Sampling**: 확률적 선택, 구현 간단

**예상 효율**:
- Active Learning 없이: 10,000-50,000 FBA 호출
- Active Learning 있이: 1,000-5,000 FBA 호출 (70-90% 감소)

---

### 2.4 모듈 D: COMETS dFBA 동적 시뮬레이션

**목적**: 미생물 군집의 시간적 동역학 시뮬레이션

**입력**:
- 각 종의 GEM 모델
- 초기 조건 (접종 비율, 배지 조성)
- 환경 파라미터 (온도, pH 등)

**출력**:
- 시간별 생체량 변화
- 시간별 대사물 농도 변화
- 군집 안정성 지표

**구현 설계**:

```python
# 의사코드
class DFBA Simulator:
    def __init__(self, comets_layout, models):
        self.layout = cometspy.layout()
        self.models = [cometspy.model(m) for m in models]
    
    def run_simulation(self, initial_biomass, max_time=168, 
                       dt=0.1, adaptive=True):
        """dFBA 시뮬레이션 실행"""
        # 초기 생체량 설정
        for model, biomass in zip(self.models, initial_biomass):
            model.initial_pop = biomass
        
        # 시뮬레이션 설정
        sim = cometspy.comets(self.layout)
        sim.parameters.maxCycles = int(max_time / dt)
        sim.parameters.timeStep = dt
        
        # 수치 안정성 설정
        if adaptive:
            # COMETS 기본: 고정 스텝
            # 커스텀: solve_ivp + BDF + 이벤트 감지
            pass
        
        sim.run()
        return self._parse_results(sim)
    
    def _parse_results(self, sim):
        """결과 파싱"""
        return {
            'biomass': sim.total_biomass,
            'media': sim.media,
            'time': sim.time_points
        }
```

**수치 안정성 대책**:

| 문제 | 원인 | 대책 |
|------|------|------|
| Stiffness | 빠른 성장 vs 느린 대사물 소비 | BDF/Radau 암시적 솔버 |
| 불연속성 | FBA basis 변화 시점 | 이벤트 감지 후 솔버 재시작 |
| Mass balance 오류 | 시간 단계가 너무 큼 | dt ≤ 0.1h, 보수적 크기 |
| 음수 농도 | Explicit 방법의 오버슈팅 | 음수 클램핑 + 암시적 방법 |

**COMETS vs 커스텀 dFBA 선택 기준**:
- **COMETS 사용**: 공간 구조 필요, 2-3종 공동 배양, 빠른 프로토타입
- **커스텀 dFBA**: 수치 정밀도 필요, 4종 이상, custom kinetics 필요

---

### 2.5 모듈 E: NSGA-II 다목적 최적화

**목적**: 미생물 군집 접종 비율 최적화

**목적 함수**:
1. **f1 (최대화)**: 타겟 대사물 생산량 (dFBA 최종 시점)
2. **f2 (최소화)**: 군집 불안정성 (생체량 변동 계수)
3. **f3 (최대화)**: 자원 효율성 (타겟/총 자원 소비)

**의사결정 변수**:
- 각 종의 접종 비율 (연속, 0.01~0.99)
- 배지 초기 농도 (선택적)

**제약 조건**:
- 접종 비율 합 = 1.0
- 최종 생체량 > 최소 임계값
- dFBA 수렴 성공

**구현**:

```python
# 의사코드 (pymoo 기반)
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

class MicrobiomeOptimization(Problem):
    def __init__(self, n_species, surrogate_model, dfba_simulator):
        super().__init__(n_var=n_species-1,  # 마지막 종 = 1 - sum(others)
                         n_obj=3,
                         n_constr=1,  # 접종 비율 합 ≤ 1
                         xl=0.01, xu=0.99)
        self.surrogate = surrogate_model
        self.dfba = dfba_simulator
    
    def _evaluate(self, X, out):
        f1, f2, f3 = [], [], []
        for x in X:
            # 접종 비율 복원
            ratios = list(x) + [1.0 - sum(x)]
            
            # 대리 모델로 빠른 예측
            pred = self.surrogate.predict(ratios)
            f1.append(pred['production'])
            f2.append(pred['instability'])
            f3.append(pred['efficiency'])
        
        out["F"] = np.column_stack([f1, f2, f3])
        out["G"] = -1.0 + np.sum(X, axis=1)  # 제약: sum ≤ 1

# 실행
algorithm = NSGA2(pop_size=100)
result = minimize(MicrobiomeOptimization(...), 
                  algorithm, 
                  termination=('n_gen', 200))
```

**계산 비용 추정**:

| 방식 | dFBA 호출 수 | 소요 시간 (8코어) |
|------|-------------|------------------|
| 대리 모델 없이 | 10,000-40,000 | 5-20시간 |
| 대리 모델 있이 | 1,000-5,000 | 30분-2시간 |
| Active Learning | 500-2,000 | 15분-1시간 |

---

### 2.6 모듈 F: TOPSIS/Pareto 의사결정 지원

**목적**: NSGA-II Pareto front에서 최종 설계점 선택

**방법론**:

**A. TOPSIS (Technique for Order of Preference by Similarity to Ideal Solution)**:
1. 정규화된 결정 행렬 구성
2. 가중 벡터 W 적용 (엔트로피 가중치 + 전문가 가중치 혼합)
3. 이상적 해(A+)와 부이상적 해(A-) 계산
4. 각 대안의 상대 근접도 계산
5. 순위 결정

**B. 가중치 민감도 분석** (필수):
- 각 가중치 +/-20% 변화에 따른 순위 안정성 검증
- 순위 역전 발생 시 Entropy weight method로 객관적 가중치 산출
- Kendall's tau로 순위 상관도 정량화

**C. 대안: Pareto front 시각화**:
- 3D scatter (f1, f2, f3)
- 2D 투영 (f1 vs f2, f1 vs f3, f2 vs f3)
- Knee point 탐지 (Pareto front의 무릎점 = 균형점)

---

## 3. 데이터 흐름도

```
┌──────────────────────────────────────────────────────────────────┐
│                                                                  │
│  ┌──────────┐  knockout   ┌──────────┐  FBA results  ┌────────┐│
│  │ Parameter│  combos     │  COBRApy │  (ground      │Training││
│  │ Space    │────────────►│  FBA     │  truth)       │ Dataset││
│  │ Sampler  │             │ (Module A)│──────────────►│  Store ││
│  └────┬─────┘             └──────────┘               └───┬────┘│
│       │                                                   │     │
│       │              ┌──────────────────────────────────┘     │
│       │              │                                        │
│       │    ┌─────────▼──────────┐                             │
│       │    │  GNN + XGBoost     │  ← Module B                │
│       │    │  Surrogate Model   │                             │
│       │    └─────────┬──────────┘                             │
│       │              │                                        │
│       │    ┌─────────▼──────────┐    uncertain     ┌────────┐ │
│       │    │  Active Learning   │    regions       │Parameter│ │
│       │    │  Loop (Module C)   │─────────────────►│Space   │ │
│       │    └─────────┬──────────┘                  │Sampler  │ │
│       │              │                              └────────┘ │
│       │    ┌─────────▼──────────┐                             │
│       │    │  COMETS dFBA       │  ← Module D                │
│       └───►│  Simulation        │                             │
│            └─────────┬──────────┘                             │
│                      │  biomass, metabolites                  │
│            ┌─────────▼──────────┐                             │
│            │  NSGA-II           │  ← Module E                │
│            │  Optimization      │                             │
│            └─────────┬──────────┘                             │
│                      │  Pareto front                         │
│            ┌─────────▼──────────┐                             │
│            │  TOPSIS/Pareto     │  ← Module F                │
│            │  Decision Support  │                             │
│            └─────────┬──────────┘                             │
│                      │  최종 설계점                            │
│                      ▼                                        │
│               실험 검증                                        │
└──────────────────────────────────────────────────────────────────┘
```

---

## 4. 개발 로드맵

### Phase 3-A: 기반 인프라 (1-2주)
| 작업 | 산출물 | 의존성 |
|------|--------|--------|
| COBRApy 병렬 FBA 파이프라인 | `fba_generator.py` | COBRApy, multiprocessing |
| GEM→이종 그래프 변환기 | `gem_to_graph.py` | COBRApy, PyTorch Geometric |
| FBA 벤치마크 수행 | `fba_benchmark_report.md` | COBRApy |

### Phase 3-B: 대리 모델 (2-3주)
| 작업 | 산출물 | 의존성 |
|------|--------|--------|
| GNN 인코더 구현 | `gnn_encoder.py` | Phase 3-A |
| XGBoost 회귀기 통합 | `surrogate_model.py` | GNN 인코더 |
| 학습/평가 파이프라인 | `train_surrogate.py` | 전체 모듈 B |

### Phase 3-C: Active Learning (1-2주)
| 작업 | 산출물 | 의존성 |
|------|--------|--------|
| Acquisition function 구현 | `acquisition.py` | 대리 모델 |
| Active Learning 루프 | `active_learning.py` | 전체 모듈 C |

### Phase 4-A: dFBA 시뮬레이션 (1-2주)
| 작업 | 산출물 | 의존성 |
|------|--------|--------|
| COMETS/cometspy 파이프라인 | `dfba_simulator.py` | cometspy |
| 수치 안정성 검증 | `stability_report.md` | dFBA 실행 결과 |
| dFBA 결과 파서 | `dfba_parser.py` | COMETS 출력 |

### Phase 4-B: 다목적 최적화 (1-2주)
| 작업 | 산출물 | 의존성 |
|------|--------|--------|
| NSGA-II 문제 정의 | `nsga2_problem.py` | pymoo + 대리 모델 |
| TOPSIS 의사결정 모듈 | `topsis_decision.py` | Pareto front |
| 통합 파이프라인 | `full_pipeline.py` | 전체 모듈 |

### Phase 4-C: 통합 검증 (1주)
| 작업 | 산출물 | 의존성 |
|------|--------|--------|
| End-to-end 테스트 | `integration_test.py` | 전체 파이프라인 |
| 결과 분석 보고서 | `final_report.md` | 실험 결과 |

**총 추정 기간**: 7-12주

---

## 5. 컴퓨팅 자원 요구사항

### Phase 3 (대리 모델)
| 자원 | 최소 | 권장 |
|------|------|------|
| GPU | RTX 3060 (12GB) | RTX 4090 (24GB) |
| CPU 코어 | 8 | 16 |
| RAM | 32GB | 64GB |
| 스토리지 | 50GB | 100GB |

### Phase 4 (dFBA + 최적화)
| 자원 | 최소 | 권장 |
|------|------|------|
| GPU | 불필요 (CPU-only 가능) | 1x GPU (대리 모델 평가용) |
| CPU 코어 | 16 | 32 |
| RAM | 64GB | 128GB |
| 스토리지 | 100GB | 200GB |

### 동시 실행 시 (Phase 3-4 통합)
| 자원 | 권장 |
|------|------|
| GPU | 2x RTX 4090 (24GB) |
| CPU | 32코어 |
| RAM | 128GB |
| 스토리지 | 300GB |

---

## 6. 리스크 매트릭스

| 리스크 | 확률 | 영향 | 완화 전략 |
|--------|------|------|-----------|
| GNN 임베딩이 동적 특성 포착 실패 | MEDIUM | HIGH | 이분 그래프 백업, GAT으로 attention 메커니즘 추가 |
| 대리 모델 일반화 실패 | MEDIUM | HIGH | Active Learning으로 적응적 데이터 수집, transfer learning |
| FLYCOP 저장소 삭제로 핵심 기능 상실 | CONFIRMED | MEDIUM | COMETS로 대체, FLYCOP 핵심 기능 재구현 |
| dFBA 수치 불안정성 | HIGH | MEDIUM | BDF 솔버, adaptive stepping, 보수적 dt |
| NSGA-II 수렴 실패 | LOW | MEDIUM | 대리 모델로 평가 비용 감소, hybrid approach |
| TOPSIS 가중치 민감도 | MEDIUM | LOW | Entropy weight, Pareto front 대안 |

---

## 7. 검증 체크리스트

### Phase 3 완료 기준
- [ ] GNN 임베딩이 노드 타입별로 의미 있는 표현 학습 (t-SNE 시각화)
- [ ] 대리 모델 R² > 0.90 (hold-out FBA 결과)
- [ ] Active Learning으로 80% 이상 FBA 호출 감소 확인
- [ ] 불확실성 추정이 실제 오차와 상관 (calibration curve)

### Phase 4 완료 기준
- [ ] COMETS dFBA 2종 공동 배양 재현 (기존 문헌 결과와 비교)
- [ ] 수치 안정성: 10회 연속 실행 모두 수렴
- [ ] NSGA-II Pareto front 도달 (100세대 이내)
- [ ] TOPSIS 순위가 전문가 판단과 일치 (Kendall's tau > 0.7)

### 통합 완료 기준
- [ ] End-to-end 파이프라인 실행 (파라미터 → 최종 설계점)
- [ ] 대리 모델 없는 NSGA-II vs 있는 NSGA-II 성능 비교
- [ ] 최종 설계점의 생물학적 타당성 검토

---

## 8. 오픈소스 도구 종합 목록

| 도구 | 버전 | 용도 | URL | 라이선스 |
|------|------|------|-----|---------|
| COBRApy | 0.31.1 | FBA 실행 | https://github.com/opencobra/cobrapy | LGPL-3.0 |
| PyTorch Geometric | latest | 이종 GNN | https://github.com/pyg-team/pytorch_geometric | MIT |
| XGBoost | latest | 회귀/분류 | https://github.com/dmlc/xgboost | Apache-2.0 |
| COMETS | 2.12.4 | dFBA 시뮬레이션 | https://github.com/segrelab/COMETS | MIT |
| cometspy | 0.6.1 | COMETS Python 인터페이스 | https://github.com/segrelab/cometspy | MIT |
| pymoo | latest | NSGA-II 최적화 | https://github.com/anyoptimization/pymoo | Apache-2.0 |
| DeepChem | latest | 분자 표현 학습 | https://github.com/deepchem/deepchem | MIT |
| scikit-learn | latest | 전처리/평가 | https://github.com/scikit-learn/scikit-learn | BSD-3 |
| scipy | latest | ODE 솔버 (BDF) | https://github.com/scipy/scipy | BSD-3 |
