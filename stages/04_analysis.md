# Stage 04 — Analysis

## Workflow
1. `docs/04_analysis.md` 가이드 확인
2. `stages/03_planning.md`에서 실험 설계 및 기준 확인
3. 데이터 전처리, 모델 학습, 평가 수행
4. 결과를 Planning 타겟과 비교
5. 산출물 → `outputs/analysis/run_XX/`에 저장
6. 아래 지식 업데이트
7. CLAUDE.md의 `current_stage` 업데이트

## Done when
- Planning 타겟 성능 달성 또는 미달 시 사용자 보고 후 방향 결정

---

## 검증된 핵심 지식

### GNN 임베딩 효과 (Exp1, run_01)
- **XGBoost-only R2=0.91 > GNN+XGBoost R2=0.82**: GNN 임베딩이 오히려 성능 저하. 137차원 knockout mask만으로 충분한 예측력 확보.
- **근거**: textbook 모델(137 genes)의 knockout mask는 이미 각 유전자의 on/off를 완전히 인코딩하므로, GNN 그래프 임베딩(32d)이 추가 정보를 제공하지 못함
- **시사점**: GNN의 가치는 knockout mask가 아닌 그래프 구조 자체의 정보(미지의 유전자 기능, pathway 간 상호작용)에서 나와야 함. 현재는 모든 유전자 기능이 알려져 있어 GNN이 redundunt

### Edge Prediction Pretraining (Exp3, run_01)
- **No pretrain R2=0.96 > Edge pretrain R2=0.91**: 사전학습이 오히려 해로움
- **근거**: edge prediction은 그래프 구조를 학습하지만, knockout→growth 예측에는 직접 관련 없음. 사전학습으로 고정된 가중치가 downstream 태스크를 방해

### Active Learning vs Random (Exp2, run_01)
- **Random R2=0.68 > AL R2=0.56**: GNN 임베딩이 무의미하면 diversity/UCB 전략도 무효
- **근거**: AL의 diversity/UCB는 임베딩 공간에서 작동하는데, 임베딩이 유용하지 않으면 전략이 random보다 나을 이유가 없음

### 양방향 엣지 필수성 (Module B 수정)
- **단방향 엣지만으로는 metabolite/gene 노드가 정보 dead-end**: HGTConv에서 dst 노드만 업데이트되므로, 모든 엣지가 reaction을 dst로 하면 metabolite/gene은 절대 업데이트 안 됨
- **해결**: reverse edges (reaction→metabolite, reaction→gene) 추가. 총 엣지 360→1036

### FBA 정답 데이터 (Module A)
- textbook 모델: single KO 137개(4.1s), double KO 500개(14.7s), random 200개(5.9s) → 총 837샘플, 24.8s
- FBA 호출당 ~30ms (GLPK solver, CPU)

---

## 다음 단계
1. ~~Step 0: 데이터 접근성 검증 (COBRApy 모델, GNN 라이브러리)~~ 완료
2. ~~Step 1-2: 기존 분석 재현 + 결과 확인~~ 완료
3. ~~Step 3: 스몰 서브셋 모듈별 구현-검증 (A→B→C)~~ 완료
4. Step 4: 설계 수정 — GNN 임베딩 전략 재검토 (ablation 필요)
5. Step 5: Full-scale Benchmark (iJO1366)

---

## Run 이력 (세부 내용은 outputs/analysis/run_XX/ 참조)
- **run_01** (2026-04-26): Module A/B/C 구현 + E2E 파이프라인. 핵심 발견: GNN 임베딩이 XGBoost-only 대비 성능 저하, pretraining 해로움, AL<Random
