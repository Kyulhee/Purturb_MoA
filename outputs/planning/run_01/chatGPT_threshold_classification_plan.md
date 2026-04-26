# 질문하는 AI에서 조건부 질문 정책으로
## LLM 기반 연구 보조 시스템의 clarification policy 관련 검토 메모

작성일: 2026-04-16  
작성 목적: 후속 테스트를 위한 중간 정리 및 연구 방향 메모

---

## 핵심 요약

[신뢰 A] 사용자의 모호한 요청을 받았을 때 무조건 답하거나 무조건 되묻는 방식보다, **기존 컨텍스트를 먼저 활용하고 필요한 경우에만 최소 질문을 하는 조건부 clarification 정책**이 더 타당하다.

[신뢰 A] 이 문제는 이미 여러 연구 축에서 다뤄지고 있다. 핵심적으로는 다음 네 영역의 교집합이다.
- asking clarification questions (ACQ)
- mixed-initiative conversational search / information seeking
- ambiguity resolution with LLMs
- preference elicitation / user modeling

[신뢰 B] 현재 문헌은 주로 일반 대화·검색·추천 환경을 다루며, **bioinformatics research assistant 맥락에서의 thresholded clarification policy**는 아직 상대적으로 비어 있다. 따라서 일반 연구축을 도메인 특화 assistant로 재구성하는 방향에 기회가 있다.

---

## 1. 문제 배경

이번 검토는 "LLM으로 바이오인포 연구를 자동화한 논문 10편을 찾아 정리해 달라"는 요청을 예시로 삼아 진행되었다. 이 사례에서 드러난 핵심 문제는 정량적 정확성보다도 **정성적 정렬 실패**였다.

구체적으로는 다음과 같은 문제가 발생할 수 있다.
- 사용자가 원하는 문헌 범위가 불명확함
- 자동화 수준이 전주기인지, 특정 단계인지 모호함
- original research / benchmark / review / tool paper가 같은 목록에 섞일 수 있음
- 사용자는 실제로는 "입문 지형 파악"을 원하는데, 시스템은 "넓은 목록 제시"로 대응할 수 있음

[신뢰 A] 이런 문제는 단순한 fact checking만으로 해결되지 않는다. 사용자의 목적, 수준, 결과 활용 맥락을 어느 정도 파악해야 한다.

---

## 2. 핵심 아이디어: 질문하는 AI가 아니라, 질문할지 결정하는 AI

[신뢰 A] 발전 방향은 "질문하는 AI" 자체보다 **언제 질문해야 하는지 판단하는 AI**로 보는 편이 더 정확하다.

권장 기본 구조는 다음과 같다.

1. **컨텍스트 기반 추론**  
   기존 대화, 사용자 프로필, 직전 작업 맥락을 바탕으로 사용자의 전문성·목적·산출물 선호를 추정한다.

2. **모호성 분해**  
   질문을 범위, 수준, 목적, 형식, 평가 기준 등으로 나누어 어디가 핵심 불확실성인지 판별한다.

3. **clarification gate**  
   불확실성이 충분히 크고, 그 불확실성이 결과 품질을 크게 바꿀 때만 질문한다.

4. **최소 질문 정책**  
   질문이 필요해도 최대 1~2개의 선택지형 질문만 한다.

5. **가정 명시 후 진행**  
   질문이 필요하지 않거나 사용자가 답하지 않으면, 합리적 기본값을 명시하고 진행한다.

[신뢰 B] 요약하면, 좋은 시스템은 always ask도 아니고 always answer도 아니다. **uncertainty-aware mixed-initiative policy**에 가깝다.

---

## 3. 왜 clarification gate가 필요한가

### 3.1 항상 바로 답하는 방식의 한계

[신뢰 A] 사용자의 질문이 겉보기에는 명확해도 실제로는 여러 해석이 가능하다. 이때 바로 답하면 그럴듯하지만 평균적인 결과에 수렴하고, 사용자가 실제로 원하는 산출물과 어긋날 가능성이 높다.

예시:
- "관련 논문 10편 찾아줘"
- "자동화 논문 정리해줘"
- "바이오인포에서 LLM 활용 연구를 보여줘"

이런 표현은 다음을 전혀 고정하지 않는다.
- 연구 논문만 볼지, 리뷰도 포함할지
- 자동화 수준을 어디까지 볼지
- 입문용 지도인지, 구현용 shortlist인지
- 최신 동향인지, 대표 고전인지

### 3.2 항상 질문하는 방식의 한계

[신뢰 A] 반대로 모든 모호한 요청에 대해 질문을 먼저 던지면 사용자가 피로해진다. 작은 형식 차이, 나중에 쉽게 수정 가능한 선택지까지 계속 되물으면 UX가 급격히 나빠진다.

### 3.3 따라서 필요한 것

[신뢰 B] 핵심은 모호성 자체가 아니라 **모호성이 결과를 얼마나 망칠 수 있는가**다. 실무적으로는 다음과 같이 볼 수 있다.

> 질문 필요도 ≈ 모호성 × 결과 영향도 × 비복구성

- **모호성**: 해석 갈래가 얼마나 많은가
- **결과 영향도**: 해석이 달라지면 결과셋이 얼마나 바뀌는가
- **비복구성**: 한 번 잘못 시작하면 뒤에서 수정 비용이 얼마나 큰가

---

## 4. 관련 연구 축

## 4.1 Asking Clarification Questions (ACQ)

[신뢰 A] clarification question 자체는 이미 독립된 연구 주제로 다뤄지고 있다. 2023년 ACL 서베이 **"A Survey on Asking Clarification Questions Datasets in Conversational Systems"**는 ACQ 데이터셋과 평가 체계를 정리하며, 비교 가능한 벤치마크 부족을 주요 문제로 지적했다.

의미:
- 네 아이디어는 완전히 새로운 발상이 아니다.
- 이미 문제 정의와 평가 체계에 대한 논의가 존재한다.

## 4.2 Mixed-initiative conversational search

[신뢰 A] 검색·정보탐색 연구에서는 clarification이 훨씬 직접적으로 다뤄진다. 사용자의 정보 요구가 불명확할 때 시스템이 먼저 clarifying question을 던지고, 이를 통해 검색 결과 품질을 높이는 방식이다.

의미:
- “논문 추천/정리” 문제는 단순 QA보다 **정보탐색**에 더 가깝다.
- 따라서 search clarification 문헌이 직접적인 참고 대상이 된다.

## 4.3 LLM 시대의 ambiguity resolution

[신뢰 A] 2024–2025년에는 **LLM이 언제 clarification을 해야 하는지**를 직접 다루는 논문들이 등장하고 있다. 대표적으로 다음 흐름이 중요하다.
- ambiguity가 있을 때 바로 대답하지 않고 interaction을 통해 해소하는 접근
- future conversation turn이나 outcome을 고려하여 질문 여부를 결정하는 접근
- uncertainty 기반으로 clarification 시점을 정하는 접근

의미:
- 네가 말한 “threshold를 넘지 못하면 먼저 질문”이라는 아이디어와 직접 연결된다.

## 4.4 Preference elicitation / user modeling

[신뢰 A] 추천 시스템과 대화형 추천 연구에서는 오래전부터 사용자 선호를 대화로 elicitation하는 문제를 다뤄 왔다. LLM 시대에는 이 흐름이 “사용자의 프로필과 목적을 먼저 추정하고, 부족할 때만 최소 질문을 하는 방식”으로 확장되고 있다.

의미:
- 네가 말한 “기존 컨텍스트에서 이 사람의 수준이나 연관성을 추론”하는 전략은 이미 연구적으로 정당화 가능한 방향이다.

---

## 5. 이 연구축을 bioinformatics assistant에 연결하면

[신뢰 B] 현재까지의 인상으로는, bioinformatics-specific clarification policy 자체를 정면으로 다루는 문헌은 아직 얇다. 대신 scientific workflow assistant, literature agent, tool-using analysis agent 쪽에서 다음 요소들이 부분적으로 나타난다.
- iterative refinement
- query decomposition
- ambiguity-aware retrieval
- conversational workflow planning

따라서 실제 연구 포지셔닝은 다음처럼 잡을 수 있다.

> 일반 LLM clarification 연구를 bioinformatics research assistant setting으로 가져와,
> “언제 질문해야 하는가”와 “무엇을 질문해야 하는가”를 도메인 특화 형태로 정의한다.

이 방향의 장점:
- 완전히 백지 연구가 아님
- 기존 일반 문헌을 baseline으로 활용 가능
- 도메인 특화 평가셋을 만들 차별점이 있음

---

## 6. 추천 문제정의

[신뢰 B] 현재 구상은 너무 넓게 잡으면 산만해질 수 있다. 초기 연구 질문은 아래처럼 좁히는 것이 좋다.

### 추천 RQ 1
**모호한 bioinformatics literature/analysis 요청에서 clarification gate가 실제로 도움이 되는가?**

비교 정책:
- always answer
- always ask
- uncertainty-gated ask

### 추천 RQ 2
**bioinformatics assistant에서 어떤 종류의 ambiguity가 결과 misalignment를 가장 크게 유발하는가?**

예시 분류:
- 문헌 유형 ambiguity
- 자동화 수준 ambiguity
- 사용자 수준 ambiguity
- 산출물 목적 ambiguity

### 추천 RQ 3
**기존 사용자 컨텍스트를 활용하면 clarification 횟수를 줄이면서도 품질을 유지할 수 있는가?**

---

## 7. 시스템 설계 방향

### 7.1 Prior-based user modeling

[신뢰 B] 질문 전에 다음을 추정한다.
- domain familiarity
- task familiarity
- evidence sensitivity
- expected output type
- likely next action

예시:
- 사용자는 바이오인포 전공자이지만, LLM-agent 문헌 지형은 처음일 수 있음
- 목적은 단순 목록이 아니라 연구 시작점 파악일 수 있음
- 이 경우 넓은 표보다 구조화된 지도가 더 유용할 가능성이 높음

### 7.2 Ambiguity decomposition

[신뢰 A] 질문을 바로 풀지 않고 먼저 아래 다섯 축으로 분해한다.
- 범위 모호성
- 수준 모호성
- 목적 모호성
- 형식 모호성
- 평가 기준 모호성

### 7.3 Clarification gate

[신뢰 B] 아래 조건 중 1개 이상이 강하게 해당되면 clarification을 권장한다.
- 문헌 유형에 따라 결과셋이 크게 달라짐
- 사용자 수준에 따라 추천 시작점이 달라짐
- high-stakes 작업인데 목적이 불명확함
- 잘못 시작하면 뒤 수정 비용이 큼

### 7.4 Minimal clarification design

[신뢰 A] 좋은 clarification은 “더 자세히 설명해 주세요”가 아니라 **선택지형 최소 질문**이어야 한다.

예시:
- "입문용 구조화 정리가 목적입니까, 구현 참고용 shortlist가 목적입니까?"
- "workflow-level 자동화 논문만 볼까요, benchmark/review도 함께 포함할까요?"

### 7.5 Safe default profile

[신뢰 B] 사용자가 답하지 않더라도 멈추지 않도록 기본 프로파일을 둔다.

예시:
- **입문 탐색형**: 넓게 보되 benchmark/review 분리
- **구현 착수형**: original research 중심, method detail 강조
- **엄밀 검토형**: inclusion/exclusion 명시, citation 기준 명시
- **실무 도입형**: 재현성·도구 의존성·운영성 강조

---

## 8. 실험 설계 초안

[신뢰 B] 바로 검증 가능한 가장 단순한 실험은 세 정책 비교다.

### 정책 A: always answer
- 질문 없이 바로 결과 생성
- 장점: 빠름
- 단점: 정성적 misalignment 위험 큼

### 정책 B: always ask
- 모호한 요청마다 먼저 질문
- 장점: 방향 정렬 강함
- 단점: 사용자 마찰 큼

### 정책 C: gated clarification
- prior 추론 후 고영향 모호성만 질문
- 최대 1~2개만 질문
- 나머지는 가정 명시 후 진행

### 평가 지표 예시
- first-response usefulness
- 후속 수정 횟수
- user correction burden
- 최종 만족도
- 답변 길이 대비 정보 밀도
- 정성적 misalignment rate

---

## 9. 이번 사례에 대한 해석

[신뢰 A] 이번 사례는 clarification gate가 실제로 작동했어야 하는 케이스로 해석된다.

이유:
- 문헌 범위 모호성이 큼
- 결과셋이 크게 달라질 수 있음
- 사용자는 해당 세부 분야에 비전문가일 수 있음
- 목적이 단순 요약이 아니라 “연구 시작점 탐색”에 가까움
- 잘못 시작하면 이후 검증 비용이 큼

이 경우의 이상적 시스템 동작은 다음과 같다.
1. 기존 맥락상 사용자가 바이오인포 연구자임을 추정
2. 다만 해당 세부 LLM 자동화 문헌에는 익숙하지 않을 수 있다고 판단
3. 목적과 범위에 대한 최소 질문 1~2개 수행
4. 답이 없으면 workflow-level original research 중심으로 기본 진행
5. benchmark/review는 별도 섹션으로 분리

---

## 10. 현재 시점의 결론

[신뢰 A] 이 주제는 이미 관련 문헌 축이 존재하는 문제이며, 연구적으로도 충분히 의미가 있다.

[신뢰 A] 다만 좋은 framing은 "질문하는 AI" 그 자체가 아니라,
**"언제 질문해야 하는지 판단하고, 질문하지 않을 때도 합리적 기본값으로 진행할 수 있는 research assistant"** 이다.

[신뢰 B] 특히 bioinformatics research assistant에 맞춘 **thresholded clarification policy**, **domain-specific ambiguity taxonomy**, **prior-based minimal questioning**은 아직 상대적으로 탐색 여지가 크다.

---

## 참고 문헌 / 참고 축

1. A Survey on Asking Clarification Questions Datasets in Conversational Systems (ACL Anthology / arXiv, 2023)  
2. Clarify When Necessary: Resolving Ambiguity Through Interaction with LMs (Findings of NAACL, 2025)  
3. Modeling Future Conversation Turns to Teach LLMs to Ask Clarifying Questions (arXiv, 2024)  
4. Asking Multimodal Clarifying Questions in Mixed-Initiative Conversational Search (arXiv, 2024)  
5. Conversational Recommender Systems (ACM Computing Surveys, 2021)  
6. Asking Clarifying Questions for Preference Elicitation with Large Language Models (Google Research)  
7. ClarQ-LLM: Evaluating Clarification and Question-Asking Capabilities of LLMs (arXiv, 2024)

---

## 검증 루틴 & 신뢰도

### Premise check
- 이번 문서는 “완전히 새로운 분야가 있는가”보다, “가까운 연구 축이 이미 존재하는가”에 초점을 맞춰 정리했다.
- 따라서 개별 논문의 세부 성과 비교보다는, 어떤 연구 전통 위에 현재 아이디어를 올릴 수 있는지를 중심으로 기술했다.

### Sanity check
- 모든 모호한 질문에 clarification을 넣으면 UX가 나빠진다.
- 반대로 아무 질문도 하지 않으면 정성적 misalignment가 반복된다.
- 따라서 조건부 clarification이라는 결론은 구조적으로 타당하다.

### Alternative path
- 대안은 질문을 직접 던지지 않고, 2~3개의 진행 옵션을 먼저 제시한 뒤 사용자가 선택하게 하는 방식이다.
- 이 접근도 clarification의 변형으로 볼 수 있으며, 일부 제품 환경에서는 더 나은 UX를 줄 수 있다.

### Uncertainty & label
- [신뢰 A]: clarification question, mixed-initiative search, preference elicitation, ambiguity resolution 관련 기존 연구가 존재한다는 점
- [신뢰 B]: bioinformatics-specific clarification policy가 아직 상대적으로 얇고 기회가 있다는 판단
- [신뢰 C]: 실제로 어떤 threshold와 정책이 최적인지는 사용자군과 도메인, 세션 길이에 따라 달라지므로 실험으로 확인해야 한다.
