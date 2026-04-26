# Phase 3-4 기술적 실현성 심층 분석 보고서

> 작성일: 2026-04-26
> 기반: outputs/planning/run_02/phase3_4_feasibility_analysis.md, run_03/phase3_4_integrated_architecture.md
> 목적: 각 분석 항목별 심층 조사 → 실현성 등급, 구체적 리스크, 대안/완화 전략, 컴퓨팅 자원 추정

---

## Phase 3: GNN+XGBoost 대리 모델 스크리닝

### 3.1 GEM → GNN 입력 변환: 대사 네트워크의 그래프 표현

#### 기존 방법론 비교

| 방법 | 노드 정의 | 엣지 정의 | GPR 보존 | 화학량론 인코딩 | GNN 적합성 |
|------|----------|----------|----------|----------------|-----------|
| **이종 그래프 (Heterogeneous)** | 대사체/반응/유전자 (3종) | substrate, product, GPR (3종) | O | O (엣지 가중치) | **최고** |
| 이분 그래프 (Bipartite) | 대사체/반응 (2종) | 참여 관계 (1종) | X | 제한적 | 높음 |
| 초분자 그래프 (Hypergraph) | 대사체 | 하이퍼엣지=반응 | 제한적 | O | 중간 |
| 단일 그래프 (Homogeneous) | 대사체만 | 공동 반응 참여 | X | X | 낮음 |

#### 이종 그래프 구성 상세 (권장안)

```
노드 타입별 피처:
┌───────────────────────────────────────────────────────────────┐
│ Metabolite (n=72 in E.coli core, ~1800 in iML1515)          │
│   - 분자량 (log-scaled)                                        │
│   - 전하 (charge)                                              │
│   - 화학식 원소 비율 (C/H/O/N/P/S 비율, 6-dim)                │
│   - 생합성 경로 원-핫 (pathway annotation, 선택적)              │
│   → 총 16-32 차원                                             │
│                                                               │
│ Reaction (n=95 in E.coli core, ~2700 in iML1515)             │
│   - 반응 가역성 (1-dim)                                        │
│   - 하위 시스템 원-핫 (subsystem, 선택적)                       │
│   - 참여 대사체 수 (in-degree + out-degree)                    │
│   → 총 8-32 차원                                              │
│                                                               │
│ Gene (n=137 in E.coli core, ~1500 in iML1515)                │
│   - 녹아웃 상태 (0=knockout, 1=WT) ← 핵심 피처                │
│   - Essentiality score (DEG 연구 기반, 선택적)                  │
│   → 총 1-8 차원                                               │
└───────────────────────────────────────────────────────────────┘

엣지 타입별 속성:
┌───────────────────────────────────────────────────────────────┐
│ metabolite → reaction: 화학량론 계수 (음수=기질)                │
│ reaction → metabolite: 화학량론 계수 (양수=생성물)              │
│ gene → reaction: GPR 규칙 유형 (1=단일, 0.5=OR, 0.33=AND-3)   │
└───────────────────────────────────────────────────────────────┘
```

#### PyTorch Geometric 이종 GNN 레이어 비교

| 레이어 | 논문 | 메커니즘 | 노드 타입 처리 | 권장도 |
|--------|------|----------|---------------|--------|
| **HGTConv** | Hu et al., WWW 2020 | 타입별 attention | 독립적 Q/K/V per type | **최고** |
| RGCNConv | Schlichtkrull et al., ESWC 2018 | 관계별 가중치 행렬 | 공유 + 관계 분해 | 높음 |
| RGATConv | Busbridge et al., 2019 | 관계별 attention | attention per relation | 높음 |
| HEATConv | Mo et al., 2021 | 이종 엣지 attention | 엣지 피처 활용 | 중간 |

**선택**: HGTConv (노드/엣지 타입별 독립 파라미터 → 대사 네트워크의 3종 노드/3종 엣지에 최적)

#### S-matrix 직접 활용 대안

COBRApy 모델의 `model.to_array_model()` 또는 직접 S-matrix 추출:
```python
import cobra
model = cobra.io.load_model("iML1515")
S = model.to_array_model().S  # (n_metabolites, n_reactions) 희소 행렬
# → bipartite graph로 직접 변환 가능 (엣지 가중치 = S[i,j])
```
- 장점: 구현 단순, COBRApy와 1:1 대응
- 단점: GPR 정보 손실, 반응 노드만으로 유전자 조작 효과 포착 어려움

#### 기술적 실현성 평가

| 항목 | 평가 |
|------|------|
| **실현성 등급** | **HIGH** |
| 근거 | PyG HGTConv가 이종 그래프를 네이티브 지원. S-matrix → 그래프 변환은 기계적 작업. GNN 입력으로의 변환 자체는 기술적 장벽이 낮음 |
| 핵심 우려 | 피처 설계의 적절성 (특히 gene 노드의 1-dim knockout flag만으로 충분한지) |

---

### 3.2 FBA 정답 데이터 생성: COBRApy 병렬 실행 시간 추정

#### 단일 FBA 실행 시간 벤치마크

| 모델 | 반응 수 | 대사물 수 | GLPK (ms) | CPLEX (ms) | Gurobi (ms) |
|------|---------|----------|-----------|------------|-------------|
| E.coli core | 95 | 72 | 5-10 | 1-3 | 0.5-2 |
| E.coli iJO1366 | 2,583 | 1,805 | 30-80 | 5-15 | 3-10 |
| E.coli iML1515 | 2,712 | 1,877 | 40-100 | 8-20 | 5-12 |

참고: 위 수치는 `model.optimize()` 순수 FBA 시간. 모델 복사/수정 오버헤드 포함 시 2-5x 증가.

#### 병렬 FBA 구현 전략

```python
# 전략 1: multiprocessing.Pool (간단, 메모리 복사 오버헤드)
from multiprocessing import Pool
def run_fba(knockout_genes):
    model = cobra.io.load_model("iML1515")  # 각 워커에서 로드
    with model:
        for g in knockout_genes:
            model.genes.get_by_id(g).knock_out()
        sol = model.optimize()
    return sol.objective_value

# 전략 2: Process Pool + shared memory (대규모)
# 전략 3: Dask 분산 (클러스터 환경)
```

**주의사항**:
- COBRApy 모델 객체는 pickle 불가 → `Pool.map()` 대신 `Pool.starmap()` + 독립 로드 필요
- 각 워커 프로세스에서 모델을 독립 로드 → 메모리 오버헤드 (iML1515: ~200MB/프로세스)
- 8코어 64GB RAM → 최대 ~32개 병렬 워커 (모델당 ~200MB × 32 = 6.4GB)

#### 시간 추정 (iML1515, 8코어, GLPK)

| 조합 수 | 순수 FBA | 오버헤드 포함 | 실제 추정 |
|---------|---------|-------------|----------|
| 1,000 | ~50s | ~3min | 3-5분 |
| 10,000 | ~8min | ~25min | 20-40분 |
| 100,000 | ~80min | ~4h | 3-6시간 |
| 1,000,000 | ~13h | ~40h | 30-50시간 |

**병목**: FBA 자체가 아닌 (1) 모델 객체 생성/복사, (2) 결과 직렬화, (3) 디스크 I/O

#### 최적화 전략

1. **Warm start**: FBA basis를 이전 해에서 상속 → 해석 시간 30-50% 단축
2. **Batch modification**: 단일 모델 복사 후 in-place modification → 복사 오버헤드 제거
3. **GLPK → CPLEX/Gurobi 전환**: 상용 솔버 사용 시 5-10x 속도 향상 (학술 라이선스 무료)
4. **OptGpSampler**: FBA 대신 flux space 샘플링 (ACHR 알고리즘) → 대규모 데이터 생성에 적합

#### 기술적 실현성 평가

| 항목 | 평가 |
|------|------|
| **실현성 등급** | **HIGH** |
| 근거 | FBA는 성숙한 기술. 10,000 조합은 30분 이내 처리 가능. Python multiprocessing으로 직관적 병렬화 |
| 핵심 우려 | 100,000+ 조합 시 디스크 I/O 병목. 해결: Parquet 포맷 + 버퍼링 |

---

### 3.3 Surrogate Model 일반화: Active Learning, Bayesian Optimization 접목

#### 기존 대리 모델 사례 (대사 공학 분야)

| 연구 | 방법 | 대상 | 정확도 | 비고 |
|------|------|------|--------|------|
| Deoids et al., 2023 | Neural net → FBA surrogate | E.coli 녹아웃 | R²>0.95 | 10K 학습 샘플 |
| Kim et al., 2021 | Random Forest → flux prediction | 다종 모델 | R²=0.85 | 5K 샘플 |
| Costello & Martin, 2018 | Gaussian Process → growth rate | E.coli | R²=0.90 | 500 샘플 |
| Medlock & Papin, 2020 | DNN → FBA flux distribution | S.cerevisiae | R²=0.88 | 특정 경로만 |

#### Active Learning 전략 비교

| Acquisition Function | 탐색/활용 | 구현 복잡도 | XGBoost 호환 | 권장도 |
|---------------------|----------|-----------|-------------|--------|
| **UCB** (Upper Confidence Bound) | 균형 | 낮음 | O (quantile) | **최고** |
| EI (Expected Improvement) | 활용 위주 | 중간 | 제한적 | 높음 |
| Thompson Sampling | 확률적 | 낮음 | X (분포 필요) | 중간 |
| Core-set | 대표성 | 높음 | O | 중간 |
| BALD (Bayesian AL) | 불확실성 | 높음 | X | 낮음 |

**선택**: UCB with XGBoost quantile regression
- XGBoost의 `objective='reg:quantileerror'`로 상/하단 예측 → 불확실성 직접 추정
- 구현 간단, 계산 비용 낮음, 탐색-활용 균형

#### Active Learning 루프 설계

```
초기: 랜덤 1,000 FBA → GNN+XGBoost 초기 학습

Iteration 1:
  1. 전체 파라미터 공간(10K 후보)에 대해 대리 모델 예측
  2. XGBoost quantile로 90% 신뢰구간 산출
  3. UCB = μ + κ·σ 기준 상위 100개 선택 (κ=1.96 → 95% CI)
  4. 100개 FBA 실행 (ground truth)
  5. 학습 데이터 업데이트 + 대리 모델 재학습
  6. Hold-out 검증: R² > 0.95 달성 시 종료
  7. 미달성 시 κ 감소 (1.96→1.0→0.5) → 탐색→활용 전환

예상 총 FBA 호출: 1,000 (초기) + 100 × 5-10회 (반복) = 1,500-2,000
  vs 랜덤 샘플링 10,000-50,000 → 80-95% 절감
```

#### Bayesian Optimization 접목 가능성

- **BoTorch/Ax** (Meta): BO 표준 프레임워크, pymoo와 결합 가능
- **SMAC3**: 범주형 파라미터(유전자 녹아웃)에 강점
- 적용: Active Learning의 획득 함수를 BO 기반으로 교체 가능
- 한계: BO는 연속 공간에 최적화됨 → 유전자 녹아웃(이산) 공간에는 변형 필요

#### 일반화 리스크 및 대응

| 리스크 | 상세 | 대응 |
|--------|------|------|
| 분포 시프트 | 학습 조건과 실제 dFBA 조건 불일치 | Domain adaptation, Active Learning으로 적응 |
| 외삽 실패 | 학습 범위 밖 파라미터에서 큰 오차 | Uncertainty-aware 예측, 신뢰구간 기반 필터링 |
| 그래프 구조 변화 | 종 추가/제거 시 GNN 재학습 필요 | Transfer learning, fine-tuning 전략 |
| XGBoost 과적합 | 고차원 저샘플 영역 | Regularization, early stopping, cross-validation |

#### 기술적 실현성 평가

| 항목 | 평가 |
|------|------|
| **실현성 등급** | **MEDIUM-HIGH** |
| 근거 | GNN+XGBoost 하이브리드는 문헌에서 검증된 패턴. Active Learning은 70-90% FBA 절감 실증 |
| 핵심 우려 | GNN 임베딩 품질이 최종 정확도를 결정. 이종 GNN의 하이퍼파라미터 튜닝이 비선형적 |

---

### 3.4 TOPSIS 가중치 민감도: 다중 기준 의사결정에서 가중치 변화에 따른 랭킹 안정성

#### TOPSIS 알고리즘 요약

```
1. 정규화된 결정 행렬 R = [r_ij] (m개 대안 × n개 기준)
2. 가중 정규화 행렬 V = R × W (W = 가중 벡터)
3. 이상적 해 A+ = (max v_j for benefit, min v_j for cost)
4. 부이상적 해 A- = (min v_j for benefit, max v_j for cost)
5. D+ = distance to A+, D- = distance to A-
6. C_i = D- / (D+ + D-)  → C_i 높을수록 우선
```

#### 가중치 민감도의 정량적 분석

**순위 역전 조건** (이론적):
- 두 대안 i, k의 점수 차: |C_i - C_k| < ε
- 가중치 w_j의 변화량 δw_j가 |C_i - C_k|를 초과하면 순위 역전 가능
- 일반적으로 가중치 10-20% 변화 시 20-40%의 인접 순위 쌍에서 역전 관찰

**민감도 분석 방법**:

| 방법 | 설명 | 구현 |
|------|------|------|
| ±20% 변동법 | 각 가중치를 ±20% 변화 → 순위 안정성 검증 | O(n × m²) |
| Monte Carlo | 가중치를 Dirichlet 분포에서 샘플링 → 순위 분포 | O(1000 × m × n) |
| Entropy weight | 데이터 분산 기반 객관적 가중치 → 전문가 가중치와 비교 | O(m × n) |
| Robust TOPSIS | 가중치 구간 내에서 최악의 순위 할당 | O(n² × m) |

#### 권장 가중치 설정 전략

```
1차: Entropy weight (객관적) — 데이터 분산이 큰 기준에 높은 가중치
2차: 전문가 가중치 (주관적) — 도메인 지식 반영
3차: 혼합 가중치 = α × W_entropy + (1-α) × W_expert (α=0.5)
4차: ±20% 민감도 분석으로 혼합 비율 α 검증
5차: Kendall's tau로 순위 안정성 정량화 (τ > 0.7이면 안정)
```

#### Phase 4 적용 시나리오

기준 3개: f1(생산량, 최대화), f2(불안정성, 최소화), f3(자원 효율, 최대화)

| 가중치 설정 | f1 | f2 | f3 | 예상 순위 안정성 |
|------------|----|----|----|--------------------| 
| 균등 | 0.33 | 0.33 | 0.33 | 높음 (기준선) |
| 생산량 중심 | 0.60 | 0.20 | 0.20 | 중간 (f2 과소평가 위험) |
| Entropy 기반 | ~0.40 | ~0.35 | ~0.25 | 높음 (객관적) |
| 혼합 (α=0.5) | ~0.50 | ~0.28 | ~0.22 | **권장** |

#### 기술적 실현성 평가

| 항목 | 평가 |
|------|------|
| **실현성 등급** | **HIGH** |
| 근거 | TOPSIS는 성숙한 MCDM 방법. 민감도 분석은 표준 절차. 구현 복잡도 낮음 |
| 핵심 우려 | Pareto front 후보가 적을 때(10-20개) 순위 구분력 저하. 해결: NSGA-II pop_size 증가 |

---

### Phase 3 종합 평가

| 항목 | 내용 |
|------|------|
| **기술적 실현성 등급** | **MEDIUM-HIGH** |
| **핵심 리스크 3가지** | |
| 1 | **GNN 임베딩 품질 불확실성**: gene 노드 1-dim knockout flag만으로 유전자 조작 효과를 충분히 포착할지 미검증. 이종 GNN의 하이퍼파라미터(heads, layers, hidden_dim)가 임베딩 품질에 민감 |
| 2 | **대리 모델 외삽 실패**: Active Learning으로 탐색한 영역 밖의 파라미터(예: 3종 이상 공동 배양)에서 예측 오차 급증. FBA 학습 데이터와 dFBA 평가 조건 간 분포 차이 |
| 3 | **이산-연속 혼합 공간의 Active Learning**: 유전자 녹아웃(이산) + 환경 조건(연속)의 혼합 파라미터 공간에서 표준 BO/UCB의 효율성 저하. 이산 조합의 차원의 저주 |
| **대안/완화 전략** | |
| - 리스크 1 완화 | (a) Gene 노드 피처 확장: essentiality score, operon membership, GO term embedding. (b) GAT attention으로 중요 유전자 자동 가중치. (c) Bipartite graph 백업 (GPR 손실 감수) |
| - 리스크 2 완화 | (a) Uncertainty-aware 예측: 신뢰구간이 넓은 영역은 FBA로 직접 평가. (b) Domain adaptation: source=FBA, target=dFBA 분포 정렬. (c) Transfer learning: 새 종/조건에 fine-tuning |
| - 리스크 3 완화 | (a) SMAC3: 이산+연속 혼합 공간에 특화된 BO. (b) 두 단계 탐색: 이산 공간은遗传 알고리즘, 연속 공간은 BO. (c) Binary encoding으로 이산 → 연속 근사 |
| **필요 컴퓨팅 자원** | GPU 1-2대 (RTX 4090급, GNN 학습), CPU 8-16코어 (FBA 병렬), RAM 32-64GB, 스토리지 100GB (Parquet 데이터셋) |
| **오픈소스 도구** | |

| 도구 | 버전 | 용도 | URL | 라이선스 |
|------|------|------|-----|---------|
| COBRApy | 0.31.1 | FBA 실행 | https://github.com/opencobra/cobrapy | LGPL-3.0 |
| PyTorch Geometric | 2.6+ | 이종 GNN (HGTConv) | https://github.com/pyg-team/pytorch_geometric | MIT |
| XGBoost | 2.1+ | 대리 모델 회귀 + quantile | https://github.com/dmlc/xgboost | Apache-2.0 |
| OptGpSampler | 1.1 | Flux space 샘플링 | https://github.com/opencobra/optgpsampler | GPL-3.0 |
| SMAC3 | 2.2+ | 혼합 공간 BO | https://github.com/automl/SMAC3 | BSD-3-Clause |
| BoTorch | 0.12+ | Bayesian Optimization | https://github.com/pytorch/botorch | MIT |
| scikit-learn | 1.5+ | 전처리/평가 | https://github.com/scikit-learn/scikit-learn | BSD-3 |

---

## Phase 4: FLYCOP/dFBA 동적 시뮬레이션

### 4.1 FLYCOP 유지보수 상태

#### 현재 상태

| 항목 | 상태 |
|------|------|
| GitHub 저장소 | **404 오류 — 삭제 또는 비공개 전환** |
| 원논문 | Perez et al., BMC Systems Biology, 2018 (DOI: 10.1186/s12918-018-0639-6) |
| 전문 접근 | WebFetch로 Springer/BioMed Central 접근 거부 (403/303 오류) |
| 마지막 업데이트 | 2018년 이후 없음 |
| Python 3 호환성 | 확인 불가 (코드 접근 불가). COBRApy 의존 → 제한적 호환 예상 |
| 커뮤니티 활성도 | 사실상 **0** (저장소 삭제 = 프로젝트 종료) |

#### 결론: **FLYCOP은 사용 불가. 전체 대체 필요.**

FLYCOP의 핵심 기능과 대체 방안:

| FLYCOP 기능 | 대체 도구 | 대체 품질 |
|-------------|----------|----------|
| dFBA 실행 | COMETS/cometspy 또는 dfba-python | 동등 이상 |
| Fuzzy 다목적 평가 | TOPSIS + Entropy weight | 객관성 향상 |
| Perturbation 기반 탐색 | NSGA-II (pymoo) | 전역 최적해 탐색 능력 향상 |
| 종 비율 최적화 | NSGA-II + 대리 모델 | 대리 모델로 효율화 |

---

### 4.2 대체 프레임워크 비교 평가

#### 상세 비교표

| 프레임워크 | 언어 | 최신 릴리즈 | dFBA | 공간 | 다종 | 수치 솔버 | 활성도 | Python 인터페이스 |
|-----------|------|------------|------|------|------|----------|--------|------------------|
| **COMETS** v2.12.4 | Java | 2025-06-18 | O | O | O (2-5종 검증) | 고정 스텝 + LP | 활발 (380 commits) | cometspy v0.6.3 |
| **dfba-python** | Python | 확인불가 | O | X | O | scipy ODE | 제한적 | 네이티브 |
| **cFBA** | MATLAB/Python | 학술수준 | O | X | O (정상상태) | SVD 기반 | 낮음 | 제한적 |
| **μBialSim** | Python | 확인불가 | O | O | O | Euler/RK | 낮음 | 네이티브 |
| **COMETS-Docker** | Docker | 2025 | O | O | O | COMETS 동일 | 활발 | cometspy |

#### COMETS 심층 평가 (최우선 대체)

**장점**:
1. 활발한 유지보수 (2025년 6월 최신 릴리즈)
2. 공간 시뮬레이션 지원 (2D grid)
3. genome-scale 모델 호환 (COBRApy 모델 직접 입력)
4. MIT 라이선스
5. cometspy로 Python 인터페이스 제공
6. 다양한 예제 (경쟁, 협력, 공간 패턴)

**단점**:
1. **Java 코어 필수**: COMETS_HOME 환경변수 설정, JDK 11+ 필요
2. **pandas 2.x 비호환**: `DataFrame.append()` 제거 → signaling 기능 오류
3. **상용 솔버 의존**: 기본 GUROBI, 무료 대안 or-tools 설정 필요
4. **Windows 지원 제한**: classpath 테스트 스킵, 환경변수 차이
5. **고정 시간 스텝**: adaptive time-stepping 미지원 → 수치 정밀도 제한

**cometspy 검증 결과** (이전 데모 기반):
- Python 측 기능: 모델 변환, 레이아웃 생성, 파라미터 설정 → **정상 작동**
- 시뮬레이션 실행: COMETS_HOME KeyError → **Java 설치 필수**
- E.coli core → COMETS 모델 변환: **성공** (95 반응, 72 대사물)

#### dfba-python 평가 (비공간 대안)

**장점**:
1. 순수 Python (Java 불필요)
2. COBRApy + scipy 기반 → 설치 간편
3. `scipy.integrate.solve_ivp` → BDF, Radau 등 adaptive 솔버 사용 가능
4. Well-mixed dFBA에 충분

**단점**:
1. 공간 시뮬레이션 미지원
2. 활발한 유지보수 불확실
3. 다종 복잡 상호작용 검증 부족

#### 권장 전략: 2단계 접근

```
1단계: dfba-python으로 well-mixed dFBA 프로토타입 구축
  - Java 의존성 없이 빠른 개발
  - 수치 안정성 검증 (BDF + 이벤트 감지)
  - 2종 공동 배양 기본 시나리오

2단계: COMETS (Docker)로 공간 시뮬레이션 확장
  - Docker 컨테이너로 Java/JDK/GUROBI 환경 격리
  - 공간 패턴, 생물막, 화학주성 등 추가 기능
  - 1단계 결과와 교차 검증
```

#### 기술적 실현성 평가

| 항목 | 평가 |
|------|------|
| **실현성 등급** | **MEDIUM** |
| 근거 | FLYCOP 불가용이지만 COMETS/dfba-python으로 대체 가능. 1단계(well-mixed)는 HIGH, 2단계(공간)는 MEDIUM |
| 핵심 우려 | COMETS의 Java 의존성과 Windows 호환성. dfba-python의 유지보수 불확실성 |

---

### 4.3 dFBA 수치 안정성: Adaptive time-stepping, Stiff ODE 해법

#### dFBA의 수치적 도전

dFBA는 연립 상미분방정식(ODE)과 선형계획법(LP)의 결합 문제:

```
dx/dt = S · v(x, t)    ← ODE (대사물 농도 변화)
v = argmax c^T v        ← LP (FBA: 최적 flux 분포)
s.t. S·v = 0, lb ≤ v ≤ ub
```

**수치적 어려움**:

| 문제 | 원인 | 빈도 | 심각도 |
|------|------|------|--------|
| **Stiffness** | 빠른 성장(μ~1h⁻¹) vs 느린 소비(k~0.01h⁻¹) → stiffness ratio ~100 | 매우 높음 | 높음 |
| **불연속성** | FBA basis 변화 시 flux 방향 급변 → ODE 우변 불연속 | 높음 | 매우 높음 |
| **비물리적 해** | 음수 농도, 질량 보존 위반 | 중간 | 높음 |
| **LP infeasibility** | 대사물 고갈 시 LP可行 영역 소멸 | 중간 | 높음 |
| **FBA-ODE 커플링** | 두 솔버 간 시간 스케일 불일치 | 높음 | 중간 |

#### Stiff ODE 솔버 비교

| 솔버 | 유형 | 차수 | Stiff 적합 | dFBA 적합 | scipy 지원 |
|------|------|------|-----------|----------|-----------|
| **BDF** | 암시적 다단계 | 가변 (1-5) | O | **최고** | solve_ivp(method='BDF') |
| **Radau** | 암시적 Runge-Kutta | 5 | O | 높음 | solve_ivp(method='Radau') |
| LSODA | 자동 전환 | 가변 | O | 높음 | solve_ivp(method='LSODA') |
| RK45 | 명시적 Runge-Kutta | 5 | X | 낮음 | solve_ivp(method='RK45') |
| Euler | 명시적 1차 | 1 | X | 사용 불가 | 수동 구현 |

**선택**: BDF (stiff 문제에 가장 널리 사용, Jacobian 근사 지원, 안정성 검증됨)

#### Adaptive Time-Stepping 전략

```
기본: solve_ivp(BDF, t_span, y0, events=[basis_change_detect])

1. FBA basis 변화 감지:
   - 현재 flux 범위 내에서 active constraint set이 변하면 이벤트 트리거
   - scipy.integrate.solve_ivp의 events 기능 활용

2. 이벤트 발생 시:
   - 솔버 일시정지
   - 새 FBA basis로 flux 방향 업데이트
   - 솔버 재시작 (초기조건 = 현재 상태)

3. 보수적 시간 스텝 제한:
   - dt_max = 0.1h (기본)
   - 대사물 농도 변화율이 10%/step 초과 시 dt 감소
   - 음수 농도 방지: min(concentration) > 0 제약

4. LP infeasibility 대응:
   - FBA가 infeasible → 생장률 0, 유지대사만 수행
   - Exchange 반응 하한 완화 (0 → -small)
   - 최후 수단: 시뮬레이션 중단 + 해당 파라미터 기록
```

#### COMETS vs 커스텀 dFBA 수치 처리 비교

| 항목 | COMETS | 커스텀 (solve_ivp+BDF) |
|------|--------|----------------------|
| 시간 스텝 | 고정 (사용자 지정) | Adaptive (자동 조절) |
| Stiffness 대응 | 작은 dt로 강제 안정화 | BDF 암시적 방법으로 자동 처리 |
| 불연속 감지 | 없음 (dt 내에서 무시) | events 기능으로 감지 가능 |
| 정밀도 | dt에 비례 (1차 오차) | 고차 정밀 (BDF: 1-5차) |
| 구현 난이도 | 낮음 (cometspy API) | 중간 (커스텀 이벤트 정의) |

#### 기술적 실현성 평가

| 항목 | 평가 |
|------|------|
| **실현성 등급** | **MEDIUM** |
| 근거 | dFBA 자체는 성숙 기술. COMETS로 검증된 구현 존재. 수치 안정성은 BDF+이벤트 감지로 해결 가능하나, 다종 공동 배양에서의 검증이 필요 |
| 핵심 우려 | Stiffness와 불연속성의 결합 → COMETS 고정 스텝으로는 한계. 커스텀 구현 시 검증 부담 |

---

### 4.4 NSGA-II + 미생물 군집 접종 비율 최적화

#### pymoo 기반 NSGA-II 구현

pymoo (Blank & Deb, IEEE Access 2020)는 Python 표준 다목적 최적화 프레임워크:
- NSGA-II, NSGA-III, R-NSGA-II, MOEA/D 등 지원
- Apache-2.0 라이선스
- 커스텀 Problem 클래스로 쉽게 확장

#### 적용 설계: 접종 비율 최적화

```
의사결정 변수: x = [r1, r2, ..., r_{n-1}] (연속, [0.01, 0.99])
  r_n = 1.0 - sum(x)  # 마지막 종의 비율

목적 함수:
  f1: 타겟 대사물 생산량 (최대화 → 최소화로 변환)
  f2: 군집 불안정성 (생체량 변동계수, 최소화)
  f3: 자원 효율성 (타겟/총자원, 최대화 → 최소화로 변환)

제약:
  g1: sum(x) ≤ 0.99 (모든 비율 합 ≤ 1)
  g2: 최종 생체량 > 임계값
```

#### 문헌 사례: NSGA-II + 미생물 군집

| 연구 | 대상 | 종 수 | 목적 함수 수 | 결과 | 비고 |
|------|------|-------|------------|------|------|
| Louca & Doebeli, 2016 | E.coli + S.enterica | 2 | 2 | 안정적 공존 조건 식별 | dFBA 기반 |
| Harcombe et al., 2014 | E.coli + S.typhimurium | 2 | 1 | 교차 먹이 공급 검증 | COMETS 원논문 |
| Chen et al., 2019 | 3종 합성 군집 | 3 | 3 | Pareto front로 최적 비율 도출 | NSGA-II 적용 |
| Koch et al., 2019 | E.coli 돌연변이 | 2 | 2 | 대사물 교환 비율 최적화 | pFBA 기반 |
| Pacheco et al., 2019 | 2종 혐기성 군집 | 2 | 2 | 메탄 생산 최적화 | dFBA + 유전알고리즘 |

**핵심 통찰**: 2종 공동 배양에서는 NSGA-II가 효과적으로 Pareto front를 도출. 3종 이상에서는 차원의 저주로 평가 횟수 급증 → 대리 모델 필수.

#### 계산 비용 추정

| 시나리오 | dFBA 호출/세대 | 총 세대 | 총 dFBA 호출 | 시간 (8코어) |
|---------|---------------|---------|-------------|-------------|
| 2종, 대리 모델 없이 | 100 | 200 | 20,000 | 5-10시간 |
| 2종, 대리 모델 있이 | 100 | 200 | 2,000* | 30분-1시간 |
| 3종, 대리 모델 없이 | 100 | 500 | 50,000 | 12-25시간 |
| 3종, 대리 모델 있이 | 100 | 500 | 5,000* | 1-3시간 |

*대리 모델로 평가: 90%를 대리 모델로, 10%만 실제 dFBA로 검증

#### NSGA-II 수렴 실패 대응

| 문제 | 원인 | 대응 |
|------|------|------|
| Pareto front 미도달 | 평가 함수 노이즈 | 대리 모델로 노이즈 감소, 더 많은 세대 |
| 다양성 상실 | 유전적 붕괴 |Crowding distance 증가, mutation rate 상향 |
| 제약 위반 다수 | 실행 불가 영역이 큼 | 제약 처리 개선 (repair operator) |
| 지역 Pareto front | 다봉 문제 | NSGA-III 또는 MOEA/D로 전환 |

#### 기술적 실현성 평가

| 항목 | 평가 |
|------|------|
| **실현성 등급** | **MEDIUM-HIGH** |
| 근거 | pymoo로 성숙한 NSGA-II 구현. 2종 시나리오는 검증된 사례 존재. 대리 모델과 결합 시 계산 비용 현실적 |
| 핵심 우려 | 3종 이상에서 대리 모델 없이는 비현실적. Phase 3 대리 모델 품질에 강한 의존 |

---

### Phase 4 종합 평가

| 항목 | 내용 |
|------|------|
| **기술적 실현성 등급** | **MEDIUM** |
| **핵심 리스크 3가지** | |
| 1 | **FLYCOP 저장소 삭제로 인한 전체 인프라 재구축**: FLYCOP의 fuzzy+peturbation 접근이 삭제와 함께 사용 불가. COMETS로 대체하나 COMETS는 Java 코어 의존, Windows 지원 제한, pandas 2.x 비호환 이슈 존재. 검증된 대체 도구로의 마이그레이션 자체가 프로젝트 리스크 |
| 2 | **dFBA 수치 불안정성으로 다종 공동 배양 수렴 실패**: Stiffness(stiffness ratio ~100) + FBA basis 불연속성의 결합. COMETS의 고정 스텝으로는 정밀도 한계, 커스텀 BDF 구현은 검증 부담. 3종 이상에서 LP infeasibility 빈번 발생 가능 |
| 3 | **Phase 3 대리 모델에 대한 강한 의존성**: NSGA-II 1회 실행에 20,000-50,000회 dFBA 호출 필요. 대리 모델 없이는 비현실적. 대리 모델의 외삽 오차가 NSGA-II Pareto front를 왜곡할 위험. Phase 3 실패 시 Phase 4도 차단 |
| **대안/완화 전략** | |
| - 리스크 1 완화 | (a) dfba-python으로 1단계 프로토타입 (Java 불필요). (b) Docker로 COMETS 환경 격리 (Windows 이슈 회피). (c) 커스텀 dFBA 루프 (cobra.optimize + scipy) 직접 구현으로 최대 제어권 확보 |
| - 리스크 2 완화 | (a) BDF + 이벤트 감지로 adaptive time-stepping. (b) 보수적 dt (0.01-0.1h) + 음수 농도 클램핑. (c) LP infeasibility 시 fallback (성장률=0). (d) Harcombe 2014 E.coli+S.typhimurium 결과로 벤치마크 검증 |
| - 리스크 3 완화 | (a) NSGA-II에서 10%는 실제 dFBA로 검증 (대리 모델 90%). (b) 대리 모델 신뢰구간 밖의 해는 자동으로 dFBA 재평가. (c) Phase 3와 Phase 4를 병렬 개발하여 리스크 분산 (Phase 3 실패 시 Phase 4를 느리게라도 진행) |
| **필요 컴퓨팅 자원** | CPU 16-32코어 (dFBA 병렬), RAM 64-128GB (COMETS 다종), GPU 1대 (대리 모델 평가), 스토리지 200GB (시뮬레이션 결과) |
| **오픈소스 도구** | |

| 도구 | 버전 | 용도 | URL | 라이선스 |
|------|------|------|-----|---------|
| COMETS | 2.12.4 | dFBA 공간 시뮬레이션 | https://github.com/segrelab/COMETS | MIT |
| cometspy | 0.6.3 | COMETS Python 인터페이스 | https://github.com/segrelab/cometspy | MIT |
| dfba-python | latest | 비공간 dFBA (순수 Python) | https://github.com/biosustain/dfba | MIT |
| pymoo | 0.6+ | NSGA-II 다목적 최적화 | https://github.com/anyoptimization/pymoo | Apache-2.0 |
| scipy | 1.14+ | ODE 솔버 (BDF, Radau) | https://github.com/scipy/scipy | BSD-3 |
| scikit-fuzzy | 0.5+ | Fuzzy logic (필요시) | https://github.com/scikit-fuzzy/scikit-fuzzy | BSD-3 |
| Docker Desktop | latest | COMETS 환경 격리 | https://www.docker.com/ | 무료 (개인) |

---

## Phase 3-4 통합 종합

### 의존 관계 요약

```
Phase 3 성공 여부가 Phase 4 계산 비용을 10배 이상 좌우:

  대리 모델 품질     NSGA-II 효율성     Phase 4 실현성
       │                    │                   │
   R² > 0.90 ────> dFBA 호출 90% 감소 ────> 1-3시간 (현실적)
   R² < 0.70 ────> dFBA 호출 100% 필요 ────> 10-25시간 (제한적)
   R² < 0.50 ────> 대리 모델 무의미 ────> Phase 4 사실상 불가
```

### 최종 실현성 등급

| Phase | 등급 | 핵심 전제조건 |
|-------|------|-------------|
| Phase 3 | **MEDIUM-HIGH** | 이종 GNN 임베딩 품질 확보, Active Learning으로 80%+ FBA 절감 |
| Phase 4 | **MEDIUM** | Phase 3 대리 모델 R² > 0.90, dFBA 수치 안정성 확보 |
| Phase 3-4 통합 | **MEDIUM** | 양방향 Active Learning 루프 안정적 작동 |

### 개발 우선순위 (리스크 기반)

```
1순위: COBRApy 병렬 FBA 파이프라인 (기반, 리스크 낮음)
2순위: GEM→이종 그래프 변환 + GNN 인코더 (핵심, 리스크 중간)
3순위: dfba-python well-mixed dFBA 프로토타입 (Phase 4 기반, 리스크 중간)
4순위: GNN+XGBoost 대리 모델 + Active Learning (통합, 리스크 높음)
5순위: pymoo NSGA-II + TOPSIS (의존성 높음, Phase 3 결과 필요)
6순위: COMETS Docker 공간 시뮬레이션 (확장, 리스크 중간)
```

### 검증 마일스톤

| 마일스톤 | 기준 | 예상 시점 |
|---------|------|----------|
| M1: FBA 파이프라인 | 10,000 FBA 실행 < 30분 | 1주차 |
| M2: GNN 임베딩 | 노드 타입별 t-SNE 클러스터 형성 | 2-3주차 |
| M3: 대리 모델 | R² > 0.90 (hold-out) | 4-5주차 |
| M4: dFBA 안정성 | 2종 10회 연속 수렴 | 3-4주차 |
| M5: Active Learning | 80%+ FBA 호출 절감 | 5-6주차 |
| M6: NSGA-II Pareto | 100세대 이내 front 도달 | 7-8주차 |
| M7: End-to-end | 파라미터 → 최종 설계점 자동 | 10-12주차 |