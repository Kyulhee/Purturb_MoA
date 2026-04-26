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

(아직 검증된 결과 없음 — 새로운 run 시작 시 이 섹션에 설계 결정의 검증 결과를 누적)

---

## 다음 단계
1. Step 0: 데이터 접근성 검증 (COBRApy 모델, GNN 라이브러리)
2. Step 1-2: 기존 분석 재현 + 결과 확인
3. Step 3: 스몰 서브셋 모듈별 구현-검증 (A→B→C→D→E→F)
4. Step 4: Baseline 비교 + 모듈 Ablation
5. Step 5: Full-scale Benchmark (iJO1366)

---

## Run 이력 (세부 내용은 outputs/analysis/run_XX/ 참조)
