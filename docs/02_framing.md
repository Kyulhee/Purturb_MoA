# Framing — Guide

## Checklist
- [ ] 문헌조사 결과(stages/01)에서 핵심 인사이트 추출
- [ ] 연구 질문이 Literature Review의 gap에서 도출되었는지 확인
- [ ] 연구 질문 명확화 (단일 문장으로 기술 가능해야 함)
- [ ] 베이스라인 수치 확정 (문헌에서 도출)
- [ ] 타겟 성능 수치 설정 (베이스라인 기반, 타당한 근거 필요)
- [ ] 사용 데이터셋 확정
- [ ] 평가 지표 확정 (AUC, F1, accuracy 등)
- [ ] 프로젝트 진행 중 연구 질문이 추가/변경된 경우, 변경 사유와 근거 기록
- [ ] 새 연구 질문 추가 시, 해당 도메인의 Literature Review를 먼저 수행 (stages/01에 gap이 있어야 질문 도출 가능)
- [ ] outputs/framing/run_XX/에 산출물 저장
- [ ] stages/02_framing.md의 Current State 업데이트

## Key Questions
- 이 질문은 어떤 기존 한계(gap)에서 비롯되었는가?
- 우리가 풀고자 하는 문제는 정확히 무엇인가?
- 성공의 기준은 무엇인가? (정량적)
- 어떤 데이터로, 어떤 지표로 측정하는가?
- 베이스라인 대비 얼마나 개선해야 의미 있는가?

## 산출물 예시
- `outputs/framing/run_01/framing_doc.md` — 문제 정의, 베이스라인, 타겟
- `outputs/framing/run_01/question_change_log.md` — 질문 추가/변경 기록 (발생 시에만)
