# Literature Review — Guide

## stages/01 서술 구조 (필수)

아래 순서로 서술. 주제별 분류가 아닌 **논리 흐름**을 따를 것:

1. **연구 배경**: 이 분야에서 풀고자 하는 근본 문제가 무엇인가?
2. **SOTA**: 현재 가장 잘 하는 방법은 무엇인가? (정량 수치 필수)
3. **베이스라인**: 재현 가능한 기준 성능은? (SOTA 아래의 실용적 기준)
4. **Gap**: SOTA/베이스라인이 해결하지 못하는 것은 무엇인가? (Framing으로 전달)

**주의사항:**
- 종료된 연구 방향은 "이전 연구 이력"으로 2-3줄 요약 후 Gap에만 인용
- 현재 연구 방향의 SOTA/베이스라인이 본문을 차지
- 서술 순서는 배경→SOTA→베이스라인→Gap의 계단 구조. 이 순서를 따르면 자연스럽게 "왜 이 연구가 필요한가"가 증명됨

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
- 이 분야에서 풀고자 하는 근본 문제는 무엇인가? (배경)
- 현재 가장 잘 하는 방법은 무엇인가? 왜 그 방법이 쓰이는가? (SOTA)
- 재현 가능한 베이스라인 성능은 어느 수준인가? (베이스라인)
- 기존 접근법의 공통된 한계(gap)는 무엇인가? — Framing의 연구 질문으로 이어짐

## 산출물 예시
- `outputs/literature_review/run_01/literature_review.md` — 전체 리뷰 문서
- `outputs/literature_review/run_01/papers/` — 개별 논문 분석
- `outputs/literature_review/run_01/gap_summary.md` — 기존 방법의 한계 요약 (Framing으로 전달)
