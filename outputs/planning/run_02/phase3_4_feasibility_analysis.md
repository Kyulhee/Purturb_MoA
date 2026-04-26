# Phase 3-4 기술적 실현성 심층 분석 보고서

> 작성일: 2026-04-26
> 분석 대상: OCT LLM XAI 프로젝트 Phase 3 (GNN+XGBoost 대리 모델 스크리닝) 및 Phase 4 (FLYCOP/dFBA 동적 시뮬레이션)

---

## Phase 3: GNN+XGBoost 대리 모델 스크리닝

### 1. GEM -> GNN 입력 변환: 대사 네트워크를 그래프로 표현

세 가지 주요 패러다임 존재:

**(A) 이종 그래프(Heterogeneous Graph) 접근 -- 권장:**
- 노드 유형: 대사체(Metabolite), 효소(Enzyme/Reaction), 유전자(Gene) 3종
- 엣지 유형: 대사체->반응(입력), 반응->대사체(출력), 유전자->반응(GPR 연관)
- PyTorch Geometric(PyG) 지원: RGCNConv, RGATConv, HGTConv, HEATConv
- 대사 네트워크의 방향성, 화학량론, GPR 규칙을 가장 충실하게 보존

**(B) 초분자 그래프(Hypergraph) 접근:**
- 반응을 하이퍼엣지로 모델링 (다중 입력->다중 출력)
- 화학량론적 계수를 엣지 가중치로 인코딩 가능

**(C) 이분 그래프(Bipartite Graph) 접근:**
- 노드: 대사체와 반응 두 종류만
- COBRApy S-matrix를 직접 그래프로 변환
- 구현이 가장 단순하나 GPR 정보 손실

### 2. FBA 정답 데이터 생성: COBRApy 병렬 실행

- COBRApy v0.31.1 (2026-03-26 최신), 활발히 유지보수 중
- COBRApy 자체 병렬 FBA 기능 없음 -> Python multiprocessing 사용
- E. coli iML1515 기준: 단일 FBA ~5-10ms (GLPK), ~1-3ms (CPLEX)
- 병목: FBA 자체보다 모델 객체 복사/수정 오버헤드

| 조합 수 | 소요 시간 (8코어, CPLEX) |
|---------|------------------------|
| 1,000   | ~1-5분                 |
| 10,000  | ~10-30분               |
| 100,000 | ~2-5시간               |

### 3. Surrogate Model 일반화: Active Learning, Bayesian Optimization 접목

- **Active Learning**: 불확실성 높은 영역 우선 탐색, 예상 FBA 호출 70-90% 감소
- **Bayesian Optimization**: Ax, BoTorch 등 활용, pymoo와 결합 가능
- **GNN+XGBoost 하이브리드**: GNN은 구조적 정보 포착, XGBoost는 테이블 데이터에 강점

### 4. TOPSIS 가중치 민감도

- 가중치 벡터의 작은 변화가 순위를 역전시킬 수 있음 (알려진 약점)
- 완화: Sensitivity analysis (+/-20%), Entropy weight method (객관적 가중치), Robust TOPSIS, Pareto front 대안

### Phase 3 종합 평가

| 항목 | 평가 |
|------|------|
| **기술적 실현성 등급** | **MEDIUM-HIGH** |
| 핵심 리스크 1 | GNN 임베딩이 동적 특성(flux 방향성, thermodynamic constraints) 포착 부족 가능 |
| 핵심 리스크 2 | 학습 데이터와 실제 생물학적 조건 간 분포 차이로 일반화 실패 |
| 핵심 리스크 3 | 이종 GNN 하이퍼파라미터 튜닝의 계산 비용 |
| **필요 컴퓨팅 자원** | GPU 1-2대 (RTX 4090급), CPU 8-16코어, RAM 32-64GB, 스토리지 100GB |

---

## Phase 4: FLYCOP/dFBA 동적 시뮬레이션

### 1. FLYCOP 유지보수 상태

- FLYCOP GitHub 저장소: **404 오류 -- 삭제되거나 비공개 전환됨**
- 2018년 이후 업데이트 없음, 사실상 **유지보수 중단/폐기** 상태
- FLYCOP 의존 연구 설계는 고위험

### 2. 대체 프레임워크 비교

| 프레임워크 | 언어 | 최신 릴리즈 | dFBA | 공간 | 활성도 | 평가 |
|-----------|------|------------|------|------|--------|------|
| **COMETS** v2.12.4 | Java+Python | 2025-06-18 | O | O | 활발 | **최우선 추천** |
| **cometspy** v0.6.1 | Python | 최신 | O | O | 활발 | COMETS 인터페이스 |
| dfba | Python | 확인불가 | O | X | 제한적 | 단순 시나리오용 |
| cFBA | Python | 확인불가 | O | X | 학술수준 | 정상상태 모델링 |
| FLYCOP | Python | **삭제됨** | O | X | **중단** | 사용 불가 |

COMETS 상세:
- Daniel Segre Lab (BU), Harcombe et al., Cell Reports, 2014
- 기능: genome-scale metabolic network + dynamic FBA + discrete diffusion
- MIT 라이선스, 380 commits
- 한계: Java 코어, cometspy 래퍼 필요

### 3. dFBA 수치 안정성

- **Stiffness**: 빠른 성장 vs 느린 대사물 소비 -> 암시적 솔버(BDF, Radau) 필요
- **Adaptive time-stepping**: 이벤트 감지로 스텝 크기 조절
- **불연속성**: FBA basis 변화 시점 -> 불연속 감지 후 솔버 재시작
- 권장: COMETS 고정 스텝 + 보수적 크기(0.01-0.1h), 커스텀은 solve_ivp+BDF+이벤트 감지

### 4. NSGA-II + 미생물 군집 접종 비율 최적화

- pymoo (Apache-2.0): NSGA-II/III, R-NSGA-II 지원
- 적용: 접종 비율 최적화 (생산량 최대화, 안정성 최소화, 자원 효율 최대화)
- 대리 모델 없이: 10,000-40,000회 dFBA 실행
- 대리 모델 활용 시: ~1,000-5,000회로 축소 가능

### Phase 4 종합 평가

| 항목 | 평가 |
|------|------|
| **기술적 실현성 등급** | **MEDIUM** |
| 핵심 리스크 1 | FLYCOP 저장소 삭제로 핵심 인프라 의존성 상실 |
| 핵심 리스크 2 | dFBA 수치 불안정성으로 다종 공동 배양 시나리오 수렴 실패 |
| 핵심 리스크 3 | NSGA-II 계산 비용 높음 (Phase 3 대리 모델 강한 의존) |
| **필요 컴퓨팅 자원** | CPU 16-32코어, RAM 64-128GB, GPU 1대, 스토리지 200GB |

---

## Phase 3-4 통합 의존성 분석

```
Phase 3 (대리 모델)          Phase 4 (동적 시뮬레이션)
        |                            |
   GNN 임베딩                  COMETS dFBA
        |                            |
   XGBoost 예측  ----적용---->  NSGA-II 평가함수
        |                            |
   Active Learning             접종 비율 최적화
        |                            |
   FBA 정답 데이터  <----요청----  미탐색 파라미터 공간
```

핵심 통찰: Phase 3과 Phase 4는 양방향 의존 관계 (Active Learning 루프)

추천 개발 순서:
1. COBRApy 병렬 FBA 파이프라인 구축
2. PyG 이종 GNN + XGBoost 대리 모델 프로토타입
3. COMETS/cometspy dFBA 파이프라인 구축
4. Active Learning 루프 구현 (Phase 3-4 연결)
5. pymoo NSGA-II + 대리 모델 다목적 최적화
6. TOPSIS/Pareto front 의사결정 지원 모듈
