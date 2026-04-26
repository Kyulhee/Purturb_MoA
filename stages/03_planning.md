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

## 연구 질문 (stages/02에서)

**"대사 네트워크의 이종 그래프에서 GNN 대리 모델이 FBA 기반 탐색 공간 축소에 기여하는가?"**

---

## 1. 파이프라인 아키텍처

```
Module A: FBA Ground Truth Generator
  → Module B: GNN+XGBoost Surrogate Model
  → Module C: Active Learning Loop (diversity → UCB)
```

의존: A→B→C 단방향. C의 AL 선택이 A의 FBA 호출을 유발하며, B가 재학습됨.

## 2. 설계 결정과 근거

| 설계 결정 | 근거 | 대안 검토 |
|-----------|------|----------|
| HGTConv 2층 이종그래프 GNN | metabolite/reaction/gene 3노드타입에 자연 대응 | GCNConv(동종만), GATConv(이종 미지원), RGCNConv(관계형, 검토 가능) |
| GNN 임베딩 + XGBoost 분리 | multi-objective loss 분리 최적화 (이전 Loss 불균형 교훈) | End-to-end GNN (loss 충돌 위험) |
| GNN fine-tuning 허용 | pretrained component fine-tuning 필수 (이전 encoder 동결 교훈) | GNN 동결 (성능 상한 제한) |
| Active Learning two-phase | R2<0에서 uncertainty 무효, diversity로 초기 탐색 | UCB from start (R2<0에서 무효) |
| COBRApy multiprocessing | FBA 병렬 실행으로 샘플 생성 가속 | 직렬 (8x 느림) |

## 3. 실험 설계

### Module A: FBA Ground Truth Generator
- **입력**: COBRApy 모델(textbook), knockout_strategy(single/double/random_subset)
- **출력**: (knockout_mask, growth_rate) 쌍
- **그래프 변환**: metabolite(72), reaction(95), gene(137) 3노드타입 이종그래프
  - 엣지: met→rxn(stoichiometry): 188+172, gene→rxn(GPR): 158
- **샘플 생성 계획**:
  - 초기: random 2,000샘플 (8-core ~24s)
  - AL 라운드당: 50샘플 (~0.6s)
  - 총 목표: 3,000-5,000샘플

### Module B: GNN+XGBoost Surrogate
- **GNN**: HGTConv 2층, hidden=32, heads=2
- **Pretraining**: autoencoder → edge prediction / contrastive loss (선택)
- **임베딩**: graph_emb(gene mean pooling) + knockout_mask concat → 169차원 (32+137)
- **XGBoost**: 임베딩 입력, objective=reg:squarederror
- **학습/평가**: 80/20 split, 5-fold CV
- **하이퍼파라미터 탐색**: GNN(hidden: [16,32,64], heads: [2,4]), XGBoost(max_depth: [3,6,9], lr: [0.01,0.1,0.3])

### Module C: Active Learning Loop
- **Phase 1 (R2 < 0.3)**: diversity 전략 — 임베딩 공간에서 최대한 먼 샘플 선택
- **Phase 2 (R2 >= 0.3)**: UCB 전략 — exploitation(high predicted growth) + exploration(high uncertainty)
- **전환 조건**: validation R2가 0.3을 안정적으로 상회 (3연속 epoch)
- **AL 라운드**: 50샘플씩, 최대 20라운드 (총 1,000 추가샘플)

## 4. 평가 설계

| 실험 | 비교 대상 | 지표 | 타겟 |
|------|----------|------|------|
| GNN 임베딩 효과 | XGBoost-only(mask) vs GNN+XGBoost | R2, RMSE | R2 > 0.5 |
| AL 탐색 효율 | Random screening vs AL | FBA 호출 수 (동일 R2 도달) | > 70% 감소 |
| AL 전환 조건 | diversity-only vs two-phase | R2 수렴 곡선 | R2=0.3에서 전환 시 수렴 가속 |
| GNN pretraining | autoencoder vs edge prediction vs contrastive | R2 | 최고 성능 pretraining 선택 |

## 5. 이전 실패 교훈의 설계 반영

| 이전 실패 | 설계 원칙 | 본 실험 반영 |
|-----------|----------|-------------|
| Loss 불균형 (50:1) | multi-objective loss는 분리 최적화 | GNN embedding loss ≠ XGBoost prediction loss |
| Encoder 동결 | pretrained component fine-tuning 허용 | GNN pretrain → fine-tune with AL |
| 평가 오류 | 평가지표는 태스크 정의와 일치 | surrogate R2 + AL 호출 감소율 (clustering metric 아님) |

## 6. 핵심 리스크와 완화

1. **샘플 부족(135→2,000+필요)** → COBRApy multiprocessing으로 2,000샘플 ~24s 생성, AL로 점진적 증가
2. **GNN pretraining 품질**(autoencoder loss 366K→365K로 거의 학습 안 됨) → edge prediction/contrastive pretraining으로 대체, ablation으로 비교
3. **AL 전환 시점 모호성** → validation R2 3연속 epoch 기준으로 자동 판단, 임계값(0.3)은 ablation으로 검증

## 7. 컴퓨팅 자원 추정

| 단계 | 예상 시간 (8-core CPU) | 비고 |
|------|----------------------|------|
| FBA 2,000샘플 생성 | ~24s | textbook 모델, multiprocessing |
| GNN pretraining | 5-15min | 모델 크기에 따라 |
| GNN+XGBoost 학습 | 1-5min | XGBoost는 빠름 |
| AL 20라운드 | 30-60min | 라운드당 재학습+FBA 50샘플 |
| 총 예상 | ~1-2h | textbook 모델 기준 |

---

## Run 이력 (세부 내용은 outputs/planning/run_XX/ 참조)
- run_01: Phase 3-4 기술실현성 심층 분석. 상세 보고서: outputs/planning/run_01/phase3_4_deep_feasibility_20260426.md
