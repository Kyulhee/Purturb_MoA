# Literature Review — Guide

## Checklist
- [ ] 검색 키워드 및 범위 정의 (질병/모델/기술 중심)
- [ ] 최소 5편 검색 (PubMed, arXiv, Semantic Scholar)
- [ ] 각 논문에서 추출: 목적, 방법, 핵심 성능 수치, 한계
- [ ] 기존 SOTA 베이스라인 수치 확보 (최소 1개)
- [ ] 각 방법론이 왜 해당 문제에서 선택되었는지(동기) 파악
- [ ] 기존 방법론의 한계(gap)를 명시 — 이 gap이 Framing으로 전달됨
- [ ] 데이터셋 가용성 확인 (공개/비공개, 접근 방법)
- [ ] outputs/literature_review/run_XX/에 산출물 저장
- [ ] stages/01_literature_review.md의 Current State 업데이트

## Key Questions
- 이 분야에서 검증된 방법론은 무엇인가? 왜 그 방법이 쓰이는가?
- 재현 가능한 베이스라인 성능은 어느 수준인가?
- 사용 가능한 공개 데이터셋은 무엇인가?
- 기존 접근법의 공통된 한계(gap)는 무엇인가? — Framing의 연구 질문으로 이어짐

## 산출물 예시
- `outputs/literature_review/run_01/literature_review.md` — 전체 리뷰 문서
- `outputs/literature_review/run_01/papers/` — 개별 논문 분석
- `outputs/literature_review/run_01/gap_summary.md` — 기존 방법의 한계 요약 (Framing으로 전달)
