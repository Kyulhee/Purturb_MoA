# Planning — Guide

## Checklist
- [ ] Framing 결과(stages/02)에서 연구 질문의 출처(gap)와 베이스라인/타겟 수치 확인
- [ ] 실험 설계 작성 (모델, 피처, 하이퍼파라미터 범위)
- [ ] 데이터 전처리 파이프라인 정의
- [ ] 교차 검증/평가 전략 수립
- [ ] 리스크 식별 및 대안 준비
- [ ] 이전 run의 실패에서 도출된 설계 원칙이 있는 경우, 이를 반영했는가?
- [ ] **사용자 컨펌 획득** (타겟 성능 + 실험 설계 승인)
- [ ] outputs/planning/run_XX/에 산출물 저장
- [ ] stages/03_planning.md의 Current State 업데이트

## stages/에 기록할 내용 범위
stages/에는 **설계 결정과 그 근거**만 기록. 상세 규칙은 CLAUDE.md stages/ 원칙 참조

## Key Questions
- 어떤 모델/방법론을 사용할 것인가? 왜?
- 데이터 분할 전략은? (train/val/test, 교차검증)
- 핵심 하이퍼파라미터와 탐색 범위는?
- 실패 시 대안은 무엇인가?
- 이 설계에서 이전 실패의 원인을 어떻게 회피하고 있는가?

## 중요
- Planning → Analysis 전환 전 **반드시 사용자 승인** 필요
- stages/03에 `✅ confirmed by user [날짜]` 기록

## 산출물 예시
- `outputs/planning/run_01/experiment_plan.md` — 실험 설계서
