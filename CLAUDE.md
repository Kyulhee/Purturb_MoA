# CLAUDE.md — Project Orchestrator

## Current Stage
`04_analysis`

## Stage Map
| Stage | Name | Guide | State File | Outputs | Object |
|-------|------|-------|------------|---------|--------|
| 01 | literature_review | docs/01_literature_review.md | stages/01_literature_review.md | outputs/literature_review/ | objects/current/idea_abstraction_card.yaml |
| 02 | framing | docs/02_framing.md | stages/02_framing.md | outputs/framing/ | objects/current/novelty_ledger.yaml |
| 03 | planning | docs/03_planning.md | stages/03_planning.md | outputs/planning/ | objects/current/experiment_contract.yaml, evaluation_validity_card.yaml |
| 04 | analysis | docs/04_analysis.md | stages/04_analysis.md | outputs/analysis/ | objects/current/result_card.yaml, validation_readiness_card.yaml, pivot_diagnosis_card.yaml |
| 05 | interpretation | docs/05_interpretation.md | stages/05_interpretation.md | outputs/interpretation/ | objects/current/claim_card.yaml |

## Rules
1. 각 단계 진입 시 해당 stages/ 문서를 먼저 읽고 워크플로우를 따를 것
2. 산출물은 outputs/ 이외에 생성하지 말 것 (단, objects/current/는 예외 — 결정 상태 추적용)
3. Planning → Analysis 전환 시 반드시 사용자 컨펌을 받을 것
4. 기존 run 결과를 덮어쓰지 말 것 (run_01, run_02... 보존)
5. Analysis 전: experiment_contract.yaml + evaluation_validity_card.yaml 확인
6. Analysis 후 (각 run): result_card.yaml + validation_readiness_card.yaml 업데이트
7. 방향 전환/루프백 전: pivot_diagnosis_card.yaml 작성 필수
8. 클레임 작성 전: claim_card.yaml + validation_readiness_card.yaml 확인
9. objects/current/ 덮어쓰기 금지 — 기존 내용이 있으면 objects/history/에 백업 후 갱신

## Loopback Rules
- Analysis → Planning: 타겟 미달 + 원인이 설계에 있을 때
- Planning → Framing: 재현 불가 + 원인이 문제 정의에 있을 때
- Framing → Literature Review: 연구 질문의 전제가 무너졌을 때
- 루프백 시: stages/는 변경 부분만 갱신, 사유를 최상단에 기록, 산출물은 새 run 번호로

## stages/ 원칙
- **정제된 합성물**: 논문 한 편의 요약이 될 정도로 최신 지식만 유지
- **이력은 outputs/에**: 취소선, 증분 패치, 변경 로그는 stages/에 두지 않음
- **Run 이력은 맨 아래 3-5줄**: 세부 내용은 outputs/ 참조
- **stages/에는 설계 결정과 그 근거만 기록**. 환경/인프라 정보, 도구 목록, 상세 코드/로그는 outputs/에 보관

## 환경
- 상세 정보: `docs/environment.md`
- System Python 3.14: 분석 스크립트 (run_13+), Ridge LOO, 통계
- ai_env (conda, Python 3.11): GEARS/CPA 학습, PyG 모델
- GPU: RTX 4060 Ti 8GB

## 실험 보고서 (사용자 요청 시)
- 가설 실패/방향 전환 시 사용자가 보고서 작성을 요청하면 `docs/07_experiment_failure_reports.md` 가이드에 따라 작성
- 저장 위치: `outputs/analysis/experiment_reports/`
- 번호 규칙: `{순번}_{키워드}.md` (한글판은 `_kr.md`)
- stages/에는 반영하지 않음

## 연구 보고서
- 가이드: `docs/08_research_report_guide.md`
- 산출물: `docs/research_report.md`
- 사용 시점: Analysis → Interpretation 전환 시, 주요 가설 검증 완료 후, 사용자 요청 시

## 작업 시작/종료 시 참조
- **Run 시작 시**: `docs/04_analysis.md` 리소스 체크 수행
- **Run 종료 시**: `docs/06_git_policy.md` 커밋/푸시 규칙 준수

## How to Resume
1. 이 파일에서 current_stage 확인
2. 해당 stages/ 문서 읽기 → Current State 파악
3. docs/ 가이드 확인
4. outputs/ 기존 산출물 파악
5. objects/current/ 현재 결정 상태 파악
6. 작업 재개
