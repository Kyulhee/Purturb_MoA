# CLAUDE.md — Project Orchestrator

## Current Stage
`05_interpretation`

## Stage Map
| Stage | Name | Guide | State File | Outputs | Object |
|-------|------|-------|------------|---------|--------|
| 01 | literature_review | docs/01_literature_review.md | stages/01_literature_review.md | outputs/literature_review/ | objects/current/idea_abstraction_card.yaml |
| 02 | framing | docs/02_framing.md | stages/02_framing.md | outputs/framing/ | objects/current/novelty_ledger.yaml |
| 03 | planning | docs/03_planning.md | stages/03_planning.md | outputs/planning/ | objects/current/experiment_contract.yaml, evaluation_validity_card.yaml, environment.yaml |
| 04 | analysis | docs/04_analysis.md | stages/04_analysis.md | outputs/analysis/ | objects/current/result_card.yaml, validation_readiness_card.yaml, pivot_diagnosis_card.yaml |
| 05 | interpretation | docs/05_interpretation.md | stages/05_interpretation.md | outputs/interpretation/ | objects/current/claim_card.yaml |

---

## Rules
1. 각 단계 진입/종료 시 해당 docs/ 가이드의 체크리스트 수행
2. 산출물은 outputs/에만 생성 (objects/current/는 예외 — 결정 상태 추적용)
3. Planning → Analysis 전환 시 반드시 사용자 컨펌
4. 기존 run 결과 덮어쓰기 금지 (run_01, run_02... 보존)
5. objects/current/ 덮어쓰기 금지 — 기존 내용 있으면 objects/history/에 백업 후 갱신
6. 루프백/방향 전환 전 pivot_diagnosis_card.yaml 작성 필수 (세부 시나리오: docs/09_loopback_protocol.md)

---

## Loopback
- 임의의 하위 단계에서 상위 단계로 루프백 가능
- AI가 도달 가능한 모든 루프백 대상을 평가하여 사유+예상 효과를 사용자에게 제시, 컨펌 받음
- 루프백 시: stages/는 변경 부분만 갱신, 사유를 최상단에 기록, 산출물은 새 run 번호로
- 세부 시나리오 + 진단 프로토콜: docs/09_loopback_protocol.md

---

## stages/ 원칙
1. **150줄 이하**: 초과 시 상세 내용을 outputs/로 이관
2. **현재 상태만**: 이전 프로젝트, 취소된 가설, 변경 전 내용은 삭제 (취소선 금지)
3. **run 이력은 최대 5줄**: 세부내용은 "outputs/{단계명}/run_XX/ 참조". run 이력과 교훈 출처는 프로젝트 전환 시에도 보존 (추적성)
4. **설계 결정과 근거만**: 환경/인프라, 도구 목록, 상세 코드/로그는 outputs/에 보관
5. **교훈 출처 분류**: same(현재 research question 하의 run) / prior(이전 프로젝트에서 전이). 판단 기준: run이 stages/02의 research question 아래에서 수행되었는가?

### 프로젝트 전환 (루프백으로 방향 변경 시)
- stages/에서 이전 프로젝트 내용 전부 삭제 (현재 research question과 직접 관련 없는 섹션). 단, run 이력과 교훈 출처는 보존
- context.yaml 초기화 + stages/ 무관 섹션 삭제
- 삭제된 내용은 outputs/에 보존되어 복구 가능
- 사용자가 이전 잔여 발견 시 AI에게 알림 → AI가 context.yaml 대조 후 삭제 제안 또는 용어 추가

---

## 용어 관리 (CONTEXT.md 패턴)
- 용어 사전: `objects/current/context.yaml` (100줄 초과 시 context_terms.yaml로 분리)
- 구조: project명 + terms(표준용어/avoid리스트) + flagged_ambiguities(모호표현→구체정의)
- 수명주기: Literature Review에서 초안 → Framing/Planning에서 추가 → Interpretation에서 구체화
- 전환 시: 기존 context를 objects/history/로 이관 후 새 프로젝트 용어로 재작성
- 충돌 시: avoid 용어 사용하면 즉시 표준 용어로 교체. 정성적 표현은 flagged_ambiguities에 판정 기준 명시

---

## Run 관리
- 모든 단계에서 run 단위 추적. run_01, run_02... 프로젝트 전체 연속 번호
- 산출물: `outputs/{단계명}/run_XX/`
- stages/에는 run 번호를 identifier로만 참조
- Run 시작 시: 해당 단계 docs/ 가이드에서 리소스 체크
- Run 종료 시: `docs/06_git_policy.md` 커밋/푸시 규칙 준수

## 환경
- 가이드: `docs/environment.md` | 실제 정보: `objects/current/environment.yaml`
- Planning 진입 시 환경 파악 후 environment.yaml 작성. AI 임의 설치 금지

## 보고서 (사용자 요청 시)
- **실험 보고서**: 가설 실패/방향 전환 시 원인 분석과 교훈 정리. 가이드: docs/07_experiment_failure_reports.md
- **연구 보고서**: 현재까지의 연구를 한눈에 파악 가능한 요약. 가이드: docs/08_research_report_guide.md

## How to Resume
1. 이 파일에서 current_stage 확인
2. objects/current/context.yaml 용어 사전 확인
3. 해당 stages/ 문서 읽기 → Current State 파악
4. docs/ 가이드 확인
5. outputs/ 기존 산출물 파악
6. objects/current/ 결정 상태 파악
7. 작업 재개
