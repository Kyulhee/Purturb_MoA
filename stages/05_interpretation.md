# Stage 05 — Interpretation

## Workflow
1. `docs/05_interpretation.md` 가이드 확인
2. `stages/04_analysis.md`에서 핵심 결과 확인
3. 결과의 도메인적 의미 해석
4. 한계점 및 후속 연구 방향 정리
5. 산출물 → `outputs/interpretation/run_XX/`에 저장
6. 아래 지식 업데이트
7. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 최종 리포트 작성 완료
- 한계점 및 후속 방향 명시됨

---

## 검증된 핵심 지식

### 이전 run 한계
- run_01 리포트는 Analysis 성능 미달 상태에서의 임시 해석
- 근본 원인(상위 단계 불충분)이 진단되지 않은 채 작성됨
- 이후 원인 분석: Framing(평가 오정의) + Planning(loss 불균형, 동결) + Literature(베이스라인 미확보)

### 원칙
- **Analysis에서 의미 있는 성능이 나온 후에만 Interpretation 진행**
- 미달 시 "Analysis 미달, 루프백"으로 처리

---

## Run 이력 (세부 내용은 outputs/interpretation/run_XX/ 참조)
- run_01: 임시 리포트 작성. 근본 원인 미진단 상태. 의미 있는 해석을 위해서는 전체 워크플로우 재진행 필요
