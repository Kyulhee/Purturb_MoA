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

**"대사 네트워크의 이종 그래프에서 GNN 대리 모델이 FBA 기반 탐색 공간 축소에 기여하는가?"**

### 도출 근거
- Literature Gap 1: 대사 네트워크 이종 그래프(metabolite/reaction/gene)에 GNN surrogate 적용 사례 부재, 현재 135샘플 R2=-0.31
- Literature Gap 2: R2<0에서 uncertainty 기반 AL이 무효, diversity→UCB 전환 조건 미검증
- Gap 1-2는 하나의 인과 체인: surrogate 예측력(R2)이 AL 유효성의 전제조건, AL이 surrogate 개선(샘플 효율)의 수단

### 질문이 답하는 것과 답하지 않는 것
- **답하는 것**: GNN 이종 그래프 임베딩이 knockout_mask 직접 사용보다 surrogate 예측력을 높이는가? AL이 random screening 대비 FBA 호출을 줄이는가?
- **답하지 않는 것**: dFBA+NSGA-II 통합 파이프라인의 수렴 보장, TOPSIS 가중치의 주관성 해소 → 후속 연구

### 평가 전략
| 평가 방식 | 지표 | 베이스라인 | 타겟 | 비고 |
|-----------|------|----------|------|------|
| Surrogate 예측력 | R2 | -0.11 (GNN-linear) | > 0.5 | FBA 정답 대비. knockout_mask만 R2=-0.31 |
| AL 탐색 효율 | FBA 호출 감소율 | 0% (random) | > 70% | 동일 R2 도달 시 호출 수 비교 |
| AL 전환 조건 | diversity→UCB 전환 R2 | 미확인 | 0.3 | UCB가 diversity를 역전하는 R2 임계값 |

### 베이스라인 계층 (surrogate 예측력)
1. Random screening (FBA 전수) — 하한
2. XGBoost-only, knockout_mask 입력 (R2=-0.31) — GNN 임베딩 없음
3. GNN-only, linear head (R2=-0.11) — XGBoost 없음
4. GNN+XGBoost, no AL — 제안 방법의 AL 제외 변형
5. GNN+XGBoost + AL (제안 방법) — 상한

### 데이터
- COBRApy textbook 모델: 95 rxn, 72 met, 137 genes (개발/검증)
- iJO1366 (E. coli): 2583 rxn, 1805 met, 1367 genes (확장 검증)
- BiGG 108개 공개 모델 (일반화 테스트)

### 이전 실패에서의 설계 원칙
| 이전 실패 | 설계 원칙 | 본 질문에서의 반영 |
|-----------|----------|-------------------|
| Loss 불균형 (50:1) | multi-objective loss는 분리 최적화 | GNN embedding loss ≠ XGBoost prediction loss |
| Encoder 동결 | pretrained component fine-tuning 필수 | GNN fine-tuning 허용 |
| 평가 오류 | 평가지표는 태스크 정의와 일치 | clustering metric → surrogate R2로 직결 |

---

## 후속 질문 (본 질문 달성 후)

**"dFBA 시뮬레이션 기반 NSGA-II 다목적 최적화가 미생물 군집 설계에서 Pareto 최적해를 안정적으로 도출하는가?"**

- 도출 근거: Literature Gap 3-4
- 전제조건: 본 질문에서 surrogate R2 > 0.5 달성 → surrogate-assisted NSGA-II 가능
- 본 질문 실패 시에도 brute-force dFBA+NSGA-II로 진행 가능하나 계산비용 급증

---

## Run 이력 (세부 내용은 outputs/framing/run_XX/ 참조)
