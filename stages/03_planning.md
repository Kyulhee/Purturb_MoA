# Stage 03 — Planning

## Workflow
1. `docs/03_planning.md` 가이드 확인
2. `stages/02_framing.md`에서 베이스라인/타겟 확인
3. 실험 설계 작성 (모델, 피처, 하이퍼파라미터, 평가 전략)
4. **사용자 컨펌 획득** (타겟 성능 + 실험 설계)
5. 산출물 → `outputs/planning/run_XX/`에 저장
6. 아래 지식 업데이트 (과거 run에서 검증된 인사이트만 통합)
7. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 실험 설계서 작성 완료
- 사용자 컨펌 획득 (`confirmed by user [날짜]`)

---

## 검증된 핵심 지식

### 1. 파이프라인 아키텍처

```
Module A: FBA Ground Truth → Module B: GNN+XGBoost Surrogate → Module C: Active Learning Loop
Module D: dFBA Simulation → Module E: NSGA-II Optimization → Module F: TOPSIS Decision
```

의존 구조: A→B→C (Phase 3), D→E→F (Phase 4), C↔D (AL↔dFBA 양방향)

### 2. 설계 결정과 근거

| 설계 결정 | 근거 | 대안 검토 |
|-----------|------|----------|
| HGTConv 2층 이종그래프 GNN | metabolite/reaction/gene 3노드타입에 자연 대응 | GCNConv(동종그래프만), GATConv(이종 미지원) |
| GNN 임베딩 + XGBoost 분리 | multi-objective loss 분리 최적화 (이전 Loss 불균형 교훈) | End-to-end GNN (loss 충돌 위험) |
| GNN fine-tuning 허용 | pretrained component fine-tuning 필수 (이전 encoder 동결 교훈) | GNN 동결 (성능 상한 제한) |
| scipy BDF/Radau dFBA 직접 구현 | FLYCOP 삭제, COMETS Java 의존 회피, 수치 안정성 확보 | COMETS(Java), Euler(5.7x 과대) |
| TOPSIS Expert(0.7/0.3) + Entropy 보조 | Expert tau=0.73, top3=100% (가장 안정) | Entropy only(tau=0.41) |
| Active Learning two-phase (diversity→UCB) | R2<0일 때 uncertainty 무효, diversity로 초기 탐색 | UCB from start (R2<0에서 무효) |
| NSGA-II 다목적 최적화 | pymoo 0.6.1.6 검증, Pareto front 30해 도출 | 단일 목적 최적화 (trade-off 불가) |

### 3. 이전 실패 교훈의 설계 반영

| 이전 실패 | 설계 원칙 | 구체적 반영 |
|-----------|----------|------------|
| Loss 불균형 (50:1) | multi-objective loss는 분리 최적화 | GNN embedding loss ≠ XGBoost prediction loss |
| Encoder 동결 | pretrained component fine-tuning 허용 | GNN autoencoder pretrain → fine-tune with AL |
| 평가 오류 | 평가지표는 태스크와 일치해야 함 | TOPSIS ranking stability (Kendall's tau) |

### 4. 컴퓨팅 자원 추정

**COBRApy FBA 병렬 실행 (8-core 기준):**
- Single knockout 137: ~2s | Double knockout 9,316: ~2min | Random 5,000: ~1min
- iJO1366 double knockout 934K: ~10h (FBA당 ~300ms 추정)

**GNN 학습:** textbook 모델 5,000샘플 → 10-30min (CPU), 5-10x 가속 (GPU)

**NSGA-II + dFBA:**
- 소규모(pop=30, gen=50): 1,500 eval × 0.78s = ~20min (textbook)
- 전체(pop=100, gen=200): 20,000 eval × 0.78s = ~4.3h (textbook), ~10-20h (iJO1366)

### 5. 핵심 리스크와 완화 전략

**Phase 3:**
1. 샘플 부족(135→500+필요) → COBRApy multiprocessing + random 5,000샘플 초기 학습 → AL로 추가
2. GNN 임베딩 품질(autoencoder 거의 학습 안 됨) → edge prediction / contrastive pretraining으로 대체
3. AL 전환 시점(UCB가 R2>0.3 이후에만 유효) → validation R2 모니터링 자동 전환

**Phase 4:**
1. NSGA-II+dFBA 계산비용(4-20h) → 소규모 초기 탐색 → surrogate-assisted NSGA-II
2. dFBA 초기 조건 민감도 → Latin Hypercube Sampling + COBRApy 기본 조건 베이스라인
3. NSGA-II 수렴 불확실성 → Hypervolume 모니터링, NSGA-III 대안 검토

### 6. 실현성 평가
- **Phase 3**: MEDIUM — 데이터 증강 + GNN pretraining 개선으로 R2 > 0.3 도달 가능, 0.5는 AL 루프에 달림
- **Phase 4**: MEDIUM — BDF/Radau 검증, TOPSIS 안정, NSGA-II+BDF 계산비용은 완화 전략 존재

---

## Run 이력 (세부 내용은 outputs/planning/run_XX/ 참조)
- run_01: Phase 3-4 기술실현성 심층 분석. 상세 보고서: outputs/planning/run_01/phase3_4_deep_feasibility_20260426.md
