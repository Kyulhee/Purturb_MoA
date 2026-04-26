# Stage 01 — Literature Review

## Workflow
1. `docs/01_literature_review.md` 가이드 확인
2. 문헌 검색 및 분석 수행
3. 산출물 → `outputs/literature_review/run_XX/`에 저장
4. 아래 지식 업데이트 (검증된 인사이트만 통합)
5. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 문헌 5편 이상 리뷰 완료
- SOTA 베이스라인 수치 최소 1개 확보
- 사용 가능한 데이터셋 확인 완료
- 기존 방법론의 한계(gap) 명시 완료

---

## 1. 연구 배경

대사 네트워크 최적화에서 설계 공간(gene knockout 조합)은 조합 폭발에 직면함. E. coli textbook 모델만 137개 유전자로 double knockout 9,316, triple knockout ~630K 조합. 실제 산업 균주(iJO1366, 1,367 genes)에서는 double knockout만 ~934K. 기존 FBA(flux balance analysis)는 각 조합을 선형계획법으로 풀어야 하므로, 전수 조사는 계산적으로 불가능.

## 2. SOTA

| 영역 | 방법 | 성능 | 비고 |
|------|------|------|------|
| 섭동 예측 (참고) | PerturbNet | R2=0.984 (sci-Plex) | GNN 기반 surrogate의 상한 참고치. 화합물-섭동 예측이며 대사 네트워크 직접 적용은 아님 |
| 대리 모델 (화학) | GNN+BO | 분자 최적화에서 5-10x 효율 | GNN 임베딩 + Bayesian Optimization. 대사 네트워크에 적용된 사례 없음 |
| dFBA 시뮬레이션 | COMETS v2.12.4 | 공간+dFBA 통합 | Java 기반, Python 래퍼(cometspy). stiff ODE 처리 검증 |
| 다목적 최적화 | NSGA-II (pymoo) | Pareto front 30해 도출 확인 | 단일 목적 대비 trade-off 가시화. 대사 네트워크 적용 사례 부재 |
| 의사결정 | TOPSIS Expert(0.7/0.3) | Kendall tau=0.73, top3=100% | 가중치 섭동 대비 가장 안정. Entropy(tau=0.41)는 보조 |

## 3. 베이스라인

| 방법 | 성능 | 재현 가능성 |
|------|------|------------|
| Random screening | FBA 전수 호출 | 기본 baseline. 100ms/call, double knockout 9,316조합 = ~15min |
| XGBoost-only (raw features) | R2=-0.31, RMSE=0.35 | knockout_mask만으로는 불가. GNN 임베딩 필요성 시사 |
| GNN-only (linear head) | R2=-0.11 | 135샘플+169차원으로 일반화 불가. 샘플 증가 시 개선 기대 |
| Euler dFBA | 5.7x 과대추정 | 사용 금소. BDF/Radau 필수 |
| COBRApy 직접 FBA | 100ms/call (textbook), 200-500ms (iJO1366) | multiprocessing으로 8-core 시 ~8x 가속 |

## 4. Gap (Framing으로 전달)

1. **GNN 대리 모델 실증 부재**: 화학/약물 분야에서 GNN surrogate는 R2>0.9 달성하나, 대사 네트워크의 이종 그래프(metabolite/reaction/gene)에 적용된 사례 없음. 현재 135샘플로 R2=-0.31, 500+샘플 + pretraining 개선 필요
2. **AL 탐색 효율 실증 부재**: Active Learning으로 FBA 호출 수를 70-90% 감소시킬 수 있다는 가설은 타당하나, R2<0에서 uncertainty 기반 AL이 무효(diversity만 유효)라는 점이 실증적으로 확인됨. R2>0.3 이후 UCB 전환 로직 미검증
3. **dFBA+다목적 최적화 통합 부재**: dFBA 수치 안정성(BDF/Radau)과 NSGA-II는 각각 검증되었으나, 두 시스템을 통합한 end-to-end 파이프라인의 계산비용(4-20h)과 수렴 보장은 미해결
4. **의사결정 객관화**: FLYCOP(저장소 삭제)의 fuzzy logic 역할을 TOPSIS+Entropy로 대체 가능하나, Expert 가중치의 주관성 문제는 잔존

---

## 이전 연구 이력
- Perturb-seq MoA: 대조 학습 기반 약물 임베딩으로 MoA 클러스터링 개선 시도. 3중 실패(Loss 불균형/encoder 동결/평가오류)로 종료. 교훈은 Gap 1-2의 설계 원칙으로 반영 (multi-objective loss 분리, fine-tuning 허용, 평가지표-태스크 일치)

---

## Run 이력 (세부 내용은 outputs/literature_review/run_XX/ 참조)
- run_01: LLM 바이오인포마틱스 10편 + Perturb-seq MoA 6모델. 갭 식별 성공, MoA 베이스라인/대조학습/데이터 미검증
- run_02: 대조학습 2편 + 섭동예측 3편 + MoA분류 3편 리뷰. 베이스라인/이론/데이터 기본 확보
