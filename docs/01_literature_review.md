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

## 검색 원칙 (Search Principles)

### 1. 문제-우선 역방향 검색 (Problem-First Reverse Search)
방법 키워드가 아닌 **문제 기술**로 검색할 것. 예: "ICM perturbation transfer"가 아니라 "predict perturbation effect in unseen cell type without target data". 같은 아이디어가 다른 용어(도메인 적응, 전이 학습, 분해 표현)로 발표될 수 있음.

**Why:** run_08에서 "ICM"/"인과 불변성"으로만 검색하여 "도메인 적응"(BuDDI), "분해 표현 전이"(XTransferCDR), "인과 최적수송"(scDRP) 용어를 사용하는 직접 경쟁자 5개를 놓침.

**How to apply:** 각 검색 세션에서 최소 3개의 다른 용어 체계로 동일 문제를 검색. 예:
- 인과추론 용어: "causal invariance", "ICM", "mechanism independence"
- 표현학습 용어: "disentangled representation", "domain adaptation", "transfer learning"
- 응용 용어: "cross-cell-type prediction", "unseen condition", "zero-shot perturbation"

### 2. 사소한 기준선 하한 (Trivial Baseline Bound)
어떤 방법을 제안하거나 문헌에서 채택하기 전에 **가장 단순한 가능한 기준선**의 성능을 먼저 추정하거나 확인할 것.

**Why:** run_08에서 Mean Shift(RPE1_ctrl_mean + K562_delta)가 R2=0.82를 달성하여, FCR+ICM R2=0.92의 실질적 우위가 0.10에 불과함이 드러남. 이 수준의 차이는 단일 세포 해상도가 아닌 평균 수준에서만 확인되면 novelty 주장이 약해짐.

**How to apply:** 문헌 리뷰 단계에서 각 문제에 대해 "가장 단순한 가능한 방법"을 명시하고, 해당 방법의 예상 성능 범위를 기록. 기존 논문에서 trivial baseline을 보고하지 않으면 직접 추정.

### 3. 개념적 등가성 감사 (Conceptual Equivalence Audit)
제안하거나 발견한 핵심 아이디어를 **다른 용어로 재표현**했을 때 기존 작업과 겹치는지 확인할 것.

**Why:** "ICM으로 z_tx를 불변하게 만든다" = "도메인 적응으로 교란 효과를 정렬한다" = "분해 표현에서 교란 잠재 변수를 도메인 불변으로 만든다" — 모두 같은 아이디어의 다른 표현. BuDDI, scDRP, XTransferCDR이 FCR-ICM과 본질적으로 같은 접근을 다른 용어로 수행함.

**How to apply:** 각 핵심 아이디어에 대해 (a) 2-3개의 다른 학술 용어로 재표현, (b) 각 용어로 다시 검색, (c) 겹치는 기존 작업이 있으면 명시. Framing 단계의 경쟁 밀도 평가로 이어짐.

## Checklist
- [ ] 검색 키워드 및 범위 정의 (질병/모델/기술 중심)
- [ ] 문제-우선 역방향 검색 수행 (최소 3개 용어 체계로 동일 문제 검색)
- [ ] 개념적 등가성 감사 수행 (핵심 아이디어를 다른 용어로 재검색)
- [ ] 최소 5편 검색 (PubMed, arXiv, Semantic Scholar)
- [ ] 각 논문에서 추출: 목적, 방법, 핵심 성능 수치, 한계
- [ ] 기존 SOTA 베이스라인 수치 확보 (최소 1개)
- [ ] 사소한 기준선 하한 추정/확인 (각 문제에 대해 가장 단순한 방법의 예상 성능)
- [ ] 각 방법론이 왜 해당 문제에서 선택되었는지(동기) 파악
- [ ] 기존 방법론의 한계(gap)를 명시 — 이 gap이 Framing으로 전달됨
- [ ] 데이터셋 가용성 확인 (공개/비공개, 접근 방법)
- [ ] outputs/literature_review/run_XX/에 산출물 저장
- [ ] stages/01_literature_review.md의 Current State 업데이트

## Key Questions
- 이 분야에서 풀고자 하는 근본 문제는 무엇인가? (배경)
- 현재 가장 잘 하는 방법은 무엇인가? 왜 그 방법이 쓰이는가? (SOTA)
- 재현 가능한 베이스라인 성능은 어느 수준인가? (베이스라인)
- **가장 단순한 가능한 방법(trivial baseline)의 성능은?** 그 위에 얼마나 개선해야 의미 있는가? (사소한 기준선 하한)
- **동일한 핵심 아이디어가 다른 용어로 이미 발표되었는가?** (개념적 등가성 감사)
- 기존 접근법의 공통된 한계(gap)는 무엇인가? — Framing의 연구 질문으로 이어짐

## 산출물 예시
- `outputs/literature_review/run_01/literature_review.md` — 전체 리뷰 문서
- `outputs/literature_review/run_01/papers/` — 개별 논문 분석
- `outputs/literature_review/run_01/gap_summary.md` — 기존 방법의 한계 요약 (Framing으로 전달)
