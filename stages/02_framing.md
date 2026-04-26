# Stage 02 — Framing

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

## 연구 질문

**"GNN+XGBoost 대리 모델 기반 Active Learning이 대사 네트워크 최적화 탐색 효율을 개선하는가?"**

도출 근거: Literature gap 3-4 — 대사 네트워크 최적화에서 FBA 기반 탐색은 조합 폭발에 직면하며, GNN 대리 모델 + Active Learning로 탐색 효율을 개선할 여지는 존재하나 실증 부재. dFBA 수치 안정성 문제와 다목적 최적화 통합 사례도 부재.

### 평가 전략
| 평가 방식 | 지표 | 비고 |
|-----------|------|------|
| Surrogate 예측력 | R2, RMSE | FBA 정답 대비 대리 모델 정확도 |
| AL 탐색 효율 | FBA 호출 수 감소율 | 동일 성능 도달 시 호출 수 비교 |
| 다목적 최적화 | Pareto front 품질 | NSGA-II 해의 수와 분포 |
| 의사결정 안정성 | TOPSIS Kendall's tau | 가중치 섭동 시 랭킹 안정성 |

### 베이스라인 계층
1. Random screening
2. GNN-only surrogate (no XGBoost)
3. XGBoost-only (no GNN, raw features)
4. GNN+XGBoost without Active Learning
5. GNN+XGBoost + Active Learning (제안 방법)

### 타겟 수치
- Surrogate R2 > 0.5 (현재 best -0.11, 개선 필요)
- FBA 호출 70-90% 감소 (vs random screening)
- Pareto front 해 30개 이상 (NSGA-II)

### 데이터
- COBRApy textbook 모델: 95 rxn, 72 met, 137 genes
- iJO1366 (E. coli): 2583 rxn, 1805 met, 1367 genes
- BiGG 108개 공개 모델

---

## 이전 연구 질문 (종료)

**"대조 학습 기반 약물 임베딩이 sci-Plex 단일 세포 섭동 데이터에서 MoA 클러스터링 품질을 개선하는가?"**

도출 근거: Literature gap 1-2. 종료 사유: Loss 불균형/encoder 동결/평가오류의 근본 원인이 설계 단계부터 존재했으며, 동일 문제에서 재시도보다 새로운 문제 설정이 합리적.

### 이전 질문에서 현재 질문으로의 교훈 이전
| 이전 실패 | 교훈 | 현재 질문 반영 |
|-----------|------|---------------|
| Loss 불균형 (50:1) | multi-objective loss는 분리 최적화 | GNN loss + XGBoost loss 분리 |
| Encoder 동결 | pretrained component fine-tuning 필수 | GNN fine-tuning 허용 |
| 평가 오류 (leave-MoA-out 분류) | 평가지표는 태스크 정의와 일치 | TOPSIS + Entropy weight 객관화 |

---

## Run 이력 (세부 내용은 outputs/framing/run_XX/ 참조)
- run_01: 이전 질문(Perturb-seq MoA) 정의. MoA 16개 확인, 평가 전략/베이스라인/타겟 설정
