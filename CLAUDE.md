# CLAUDE.md — Project Orchestrator

## Current Stage
`03_planning`

## Stage Map
| Stage | Name | Guide | State File | Outputs |
|-------|------|-------|------------|---------|
| 01 | literature_review | docs/01_literature_review.md | stages/01_literature_review.md | outputs/literature_review/ |
| 02 | framing | docs/02_framing.md | stages/02_framing.md | outputs/framing/ |
| 03 | planning | docs/03_planning.md | stages/03_planning.md | outputs/planning/ |
| 04 | analysis | docs/04_analysis.md | stages/04_analysis.md | outputs/analysis/ |
| 05 | interpretation | docs/05_interpretation.md | stages/05_interpretation.md | outputs/interpretation/ |

## Rules
1. 각 단계 진입 시 해당 stages/ 문서를 먼저 읽고 워크플로우를 따를 것
2. 산출물은 outputs/ 이외에 생성하지 말 것
3. Planning → Analysis 전환 시 반드시 사용자 컨펌을 받을 것
4. 기존 run 결과를 덮어쓰지 말 것 (run_01, run_02... 보존)

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

## 작업 시작/종료 시 참조
- **Run 시작 시**: `docs/04_analysis.md` 리소스 체크 수행
- **Run 종료 시**: `docs/06_git_policy.md` 커밋/푸시 규칙 준수

## How to Resume
1. 이 파일에서 current_stage 확인
2. 해당 stages/ 문서 읽기 → Current State 파악
3. docs/ 가이드 확인
4. outputs/ 기존 산출물 파악
5. 작업 재개
