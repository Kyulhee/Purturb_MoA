# LLM으로 바이오인포매틱스 연구를 자동화한 논문 정리 (10편)

> 검색일: 2026-04-16 | 검색소스: Semantic Scholar, PubMed, arXiv

---

## 1. AutoBA: An AI Agent for Fully Automated Multi-Omic Analyses

| 항목 | 내용 |
|------|------|
| **저자** | Zhou Juexiao, Zhang Bin, Li Guowei, Chen Xiuying 외 |
| **연도/저널** | 2024, *Advanced Science* |
| **DOI** | 10.1002/advs.202407094 |
| **인용수** | 48 |

- **목적**: 멀티오믹스 데이터 분석을 완전 자동화하는 LLM 기반 AI 에이전트 개발
- **방법**: GPT-4 등 LLM 백엔드 기반 자율 에이전트. 최소 사용자 입력으로 단계별 분석 계획 수립 및 실행. Automated Code Repair(ACR) 메커니즘으로 코드 오류 자동 수정. 온라인/로컬 LLM 백엔드 지원으로 데이터 보안 확보
- **핵심 성과**: 사전 정의된 파이프라인이 아닌 입력 데이터에 따라 분석 과정을 자율 설계. ACR 메커니즘으로 ChatGPT 및 오픈소스 LLM 대비 end-to-end 안정성 향상. 기존 온라인 바이오인포 서비스 대비 유연성 및 적응성 우수
- **의의**: 바이오인포매틱스 분석의 "코딩 없이 자연어만으로" 접근 가능하게 만든 선도적 시스템
- **한계**: 복잡한 de novo 분석보다는 기존 도구 조합에 의존; 대규모 데이터셋에서 실행 시간 제약

---

## 2. DrBioRight 2.0: An LLM-Powered Bioinformatics Chatbot for Cancer Functional Proteomics

| 항목 | 내용 |
|------|------|
| **저자** | Liu Wei, Li Jun, Tang Yitao, Zhao Yining 외 |
| **연도/저널** | 2025, *Nature Communications* |
| **DOI** | 10.1038/s41467-025-57430-4 |
| **인용수** | 33 |

- **목적**: 대규모 암 기능 프로테오믹스 분석을 자동화하는 LLM 기반 챗봇 시스템 개발
- **방법**: RPPA(Reverse Phase Protein Array) 기반 ~8,000명 환자 샘플의 프로테오믹스 데이터 활용. LLM 기반 자연어 인터페이스로 복잡한 분석 쿼리 처리
- **핵심 성과**: 암 바이오마커 및 치료 타겟 발견을 위한 포괄적 프로테오믹스 리소스 구축. 비전문가도 자연어로 프로테오믹스 데이터 분석 가능
- **의의**: 암 연구에서 LLM이 실제 대규모 임상 데이터 분석을 자동화한 대표적 사례
- **한계**: 프로테오믹스 특화 도메인으로 확장성 제한; 환각 가능성 완전 배제 어려움

---

## 3. BixBench: A Comprehensive Benchmark for LLM-based Agents in Computational Biology

| 항목 | 내용 |
|------|------|
| **저자** | Ludovico Mitchener, Jon M. Laurent, Alex Andonian, Benjamin Tenmann 외 |
| **연도/저널** | 2025, *arXiv preprint* |
| **arXiv** | 2503.00096 |
| **인용수** | 49 |

- **목적**: LLM 에이전트의 계산생물학 실무 능력을 측정하는 종합 벤치마크 구축
- **방법**: 50개 이상의 실제 생물학적 데이터 분석 시나리오 + ~300개 서술형 문제. GPT-4o, Claude 3.5 Sonnet 평가. 오픈소스 에이전트 프레임워크 공개
- **핵심 성과**: 최첨단 프론티어 모델도 서술형 정답률 **17%**에 불과; 객관식에서는 **랜덤 수준**. 긴 다단계 분석 궤적 수행 및 미묘한 결과 해석 능력이 현저히 부족함을 입증
- **의의**: LLM 기반 바이오인포 자동화의 현재 한계를 가장 명확히 보여주는 벤치마크
- **한계**: 평가 기준이 서술형 정답 매칭에 의존; 새로운 모델에 대한 지속적 업데이트 필요

---

## 4. BioResearcher: From Intention to Implementation — Automating Biomedical Research via LLMs

| 항목 | 내용 |
|------|------|
| **저자** | Yi Luo, Linghang Shi, Yihao Li, Aobo Zhuang, Yeyun Gong 외 |
| **연도/저널** | 2024, *Science China Information Sciences* |
| **DOI** | 10.1007/s11432-024-4485-0 / arXiv: 2412.09429 |
| **인용수** | 30 |

- **목적**: 바이오메디컬 dry lab 연구 전 과정을 자동화하는 end-to-end 시스템 개발
- **방법**: 다중 에이전트 아키텍처 — 검색 에이전트, 문헌 처리 에이전트, 실험 설계 에이전트, 프로그래밍 에이전트. 계층적 학습(heierarchical learning)으로 복잡 작업 분해. LLM 기반 리뷰어로 in-process 품질 관리
- **핵심 성과**: 8개 미해결 연구 목표에서 평균 실행 성공률 **63.07%**. 생성된 프로토콜이 일반 에이전트 시스템 대비 5개 품질 지표에서 평균 **22.0%** 우수
- **의의**: 문헌 검색→실험 설계→코드 구현→결과 검증 전 사이클 자동화 시도. 바이오메디컬 연구 자동화의 가장 포괄적 엔드투엔드 시스템
- **한계**: wet lab 실험 불가; 37% 실패율; 복잡한 다학제 요구사항에서 논리적 오류 발생 가능

---

## 5. DrugAgent: Automating AI-aided Drug Discovery Programming through LLM Multi-Agent Collaboration

| 항목 | 내용 |
|------|------|
| **저자** | Sizhe Liu, Yizhou Lu, Siyu Chen, Xiyang Hu, Jieyu Zhao 외 |
| **연도/저널** | 2024, *arXiv preprint* |
| **arXiv** | 2411.15692 |
| **인용수** | 36 |

- **목적**: 신약 발견 ML 프로그래밍을 자동화하는 다중 에이전트 프레임워크 개발
- **방법**: LLM Planner(고수준 아이디어 수립) + LLM Instructor(도메인 지식 통합 및 구현) 협업. 3가지 대표적 신약 발견 태스크로 케이스 스터디
- **핵심 성과**: Drug-Target Interaction(DTI) 태스크에서 ReAct 대비 ROC-AUC **4.92%** 상대 개선. 기존 베이스라인 대비 일관된 성능 향상
- **의의**: 신약 발견이라는 고도로 전문적인 도메인에서 LLM 에이전트가 실제로 프로그래밍을 자동화한 선도 사례
- **한계**: 3개 태스크로 제한된 평가; 복잡한 분자 생성/최적화 태스크에서는 미검증; 도메인 지식 누락 시 오류 발생

---

## 6. PROTEUS: Automating Exploratory Proteomics Research via Language Models

| 항목 | 내용 |
|------|------|
| **저자** | Ning Ding, Shang Qu, Linhai Xie, Yifei Li, Zaoqu Liu 외 |
| **연도/저널** | 2024, *arXiv preprint* |
| **arXiv** | 2411.03743 |
| **인용수** | 9 |

- **목적**: 원시 프로테오믹스 데이터에서 과학적 발견을 완전 자동화하는 시스템 개발
- **방법**: LLM 기반 계층적 계획 수립 → 전문 바이오인포 도구 실행 → 반복적 워크플로우 정제. 12개 프로테오믹스 데이터셋(면역세포, 종양, single-cell/bulk) 평가
- **핵심 성과**: **191개** 과학적 가설 자동 생성. 5개 지표 LLM 자동 평가 + 인간 전문가 상세 리뷰 모두에서 일관되게 신뢰할 수 있는 결과. 기존 문헌과 잘 부합하면서도 새로운 평가 가능한 가설 제안
- **의의**: 원시 데이터 → 연구 목표 설정 → 분석 → 가설 생성 전 파이프라인을 인간 개입 없이 수행한 독보적 시스템
- **한계**: 12개 데이터셋으로 제한된 검증; 생성된 가설의 실험적 검증은 별도 수행 필요

---

## 7. Agentic Bioinformatics: Streamline Automated Biomedical Discoveries

| 항목 | 내용 |
|------|------|
| **저자** | Zhou Juexiao, Jiang Jindong, Han Zhongyi, Wang Zijian, Gao Xin |
| **연도/저널** | 2025, *Briefings in Bioinformatics* |
| **DOI** | 10.1093/bib/bbaf505 |
| **인용수** | 5 |

- **목적**: "에이전트 바이오인포매틱스(Agentic Bioinformatics)"라는 새 연구 분야의 원리, 방법론, 응용 종합 정리
- **방법**: 리뷰/비전 논문. LLM 기반 자율 에이전트의 핵심 원칙, 진화하는 방법론, 다양한 응용 분야 체계적 분석
- **핵심 내용**: 개인 맞춤 의학, 신약 발견, 합성생물학 등 핵심 영역에서 에이전트 프레임워크 통합 분석. 자율적이고 적응적인 에이전트가 생물학적 데이터셋의 자기 주도 탐색을 가능하게 하는 방식 설명. 윤리적, 기술적, 확장성 과제 식별
- **의의**: LLM 에이전트 기반 바이오인포매틱스를 하나의 독립된 연구 분야로 정의한 선도적 비전 논문
- **한계**: 리뷰 논문으로 실험적 검증 없음; 구체적 벤치마크나 정량적 비교 부족

---

## 8. GPTCelltype: Assessing GPT-4 for Cell Type Annotation in Single-Cell RNA-seq

| 항목 | 내용 |
|------|------|
| **저자** | Hou Wenpin, Ji Zhicheng |
| **연도/저널** | 2024, *Nature Methods* |
| **DOI** | 10.1038/s41592-024-02235-4 |
| **인용수** | 199 |

- **목적**: GPT-4를 활용한 single-cell RNA-seq 세포타입 주석 달기 자동화 평가
- **방법**: 마커 유전자 정보를 입력으로 GPT-4가 세포타입 주석 생성. 수백 개의 조직 및 세포타입에서 수동 주석과 비교 평가. R 패키지 GPTCelltype 개발
- **핵심 성과**: GPT-4 생성 주석이 수동 주석과 **강한 일치도** 보임. 세포타입 주석에 필요한 노력과 전문성을 상당히 줄일 수 있음을 입증
- **의의**: LLM이 바이오인포매틱스 핵심 태스크(세포타입 주석)에서 실용적 도구로 기능할 수 있음을 *Nature Methods*에 입증한 가장 인용된 논문 (199회)
- **한계**: 마커 유전자 선택에 의존; 희귀 세포타입에서 정확도 저하 가능; GPT-4 API 비용

---

## 9. Bioinfo-Bench: A Simple Benchmark Framework for LLM Bioinformatics Skills Evaluation

| 항목 | 내용 |
|------|------|
| **저자** | Qiyuan Chen, Cheng Deng |
| **연도/저널** | 2023, *bioRxiv preprint* |
| **DOI** | 10.1101/2023.10.18.563023 |
| **인용수** | 11 |

- **목적**: LLM의 바이오인포매틱스 지식 및 데이터 마이닝 능력 평가 벤치마크 구축
- **방법**: 지식 습득(knowledge acquisition), 지식 분석(knowledge analysis), 지식 응용(knowledge application) 3관점에서 체계적 데이터 수집. ChatGPT, LLaMA, Galactica 평가
- **핵심 성과**: LLM은 **지식 습득**(기억/재생)에 뛰어나지만, **실무 문제 해결** 및 **뉘앙스 있는 지식 추론** 능력은 제한적임을 입증
- **의의**: LLM의 바이오인포매틱스 능력을 체계적으로 평가한 최초의 벤치마크 중 하나. "챗봇 이상의 역할" 기대와 실제 능력 간 격차를 정량화
- **한계**: 프로젝트 진행 중으로 공개 자료 제한; 실제 코드 실행 능력 평가 미포함

---

## 10. ICGI: Cancer Gene Identification via Causal Prompting LLM + Omics Causal Inference

| 항목 | 내용 |
|------|------|
| **저자** | Zeng Haolong, Yin Chaoyi, Chai Chunyang, Wang Yuezhu, Dai Qi, Sun Huiyan |
| **연도/저널** | 2025, *Briefings in Bioinformatics* |
| **DOI** | 10.1093/bib/bbaf113 |
| **인용수** | 5 |

- **목적**: 인과적 프롬프팅 LLM + 오믹스 데이터 기반 인과 추론을 결합한 암 유전자 식별 프레임워크 개발
- **방법**: LLM에 인과성 문맥 단서(causality contextual cues)를 프롬프트로 제공 + 데이터 기반 인과적 특징 선택 모듈 결합. TCGA 6개 암 종류의 전사체 데이터 평가
- **핵심 성과**: 기존 SOTA 방법 대비 암 유전자 식별에서 **우수한 성능** (암/정상 샘플 구분). 온라인 서비스 플랫폼 구축 — 유전자 + 암 종류 입력 시 자동으로 인과적 역할 분석 결과 제공
- **의의**: LLM의 인과적 추론 능력을 오믹스 분석에 활용한 혁신적 접근. 프롬프트 엔지니어링이 바이오인포매틱스 분석 품질에 미치는 영향 입증
- **한계**: 현재 LLM은 모든 오믹스 수준의 포괄적 정보를 포착하지 못함; 6개 암 종류로 제한된 검증

---

## 종합 비교표

| # | 시스템 | 연도 | 저널 | 인용 | 자동화 범위 | 에이전트 유형 | 평가 방식 |
|---|--------|------|------|------|-------------|--------------|-----------|
| 1 | AutoBA | 2024 | Adv. Sci. | 48 | 멀티오믹스 분석 | 단일 에이전트 + ACR | 정성/정량 |
| 2 | DrBioRight 2.0 | 2025 | Nat. Commun. | 33 | 암 프로테오믹스 | 챗봇 인터페이스 | 정량 |
| 3 | BixBench | 2025 | arXiv | 49 | 계산생물학 전반 | 벤치마크 | 300문제 서술형 |
| 4 | BioResearcher | 2024 | Sci. China | 30 | 바이오메디컬 연구 전체 | 다중 에이전트 | 63.07% 성공률 |
| 5 | DrugAgent | 2024 | arXiv | 36 | 신약 발견 ML | 다중 에이전트 | ROC-AUC 비교 |
| 6 | PROTEUS | 2024 | arXiv | 9 | 프로테오믹스 발견 | 계층적 에이전트 | 191 가설 평가 |
| 7 | Agentic Bioinfo | 2025 | Brief. Bioinf. | 5 | 리뷰/비전 | — | — |
| 8 | GPTCelltype | 2024 | Nat. Methods | 199 | 세포타입 주석 | LLM 직접 호출 | 수동 주석 비교 |
| 9 | Bioinfo-Bench | 2023 | bioRxiv | 11 | LLM 역량 평가 | 벤치마크 | 3관점 평가 |
| 10 | ICGI | 2025 | Brief. Bioinf. | 5 | 암 유전자 식별 | LLM + 인과추론 | SOTA 비교 |

---

## 주요 시사점

1. **다중 에이전트가 대세**: AutoBA, BioResearcher, DrugAgent, PROTEUS 모두 LLM 에이전트(또는 다중 에이전트) 아키텍처 채택. 단일 LLM 호출보다 도구 사용+계획 수립+반복 정제 파이프라인이 효과적
2. **아직 초기 단계**: BixBench(17% 정답률)와 Bioinfo-Bench(실무 추론 한계)가 보여주듯, 완전 자율 바이오인포 분석은 요원. 현재는 "인간-인-더-루프" 보조 도구로 활용이 현실적
3. **가장 성공적 적용**: GPTCelltype(199 인용) — 특정 태스크(세포타입 주석)에 LLM을 정밀 적용한 사례가 가장 실용적 성과
4. **평가의 어려움**: 생성된 가설/프로토콜의 품질 평가가 주요 과제. 자동 평가(LLM-as-judge)와 인간 전문가 평가의 일치도 검증 필요
5. **도메인 특화 vs 범용**: 신약 발견(DrugAgent), 프로테오믹스(PROTEUS), 프로테오믹스 챗봇(DrBioRight) 등 도메인 특화 시스템이 범용 시스템보다 성능 우수

---

*이 리포트는 Semantic Scholar, PubMed, arXiv 검색 결과와 메타데이터/초록을 기반으로 작성되었습니다. PDF 전문 파싱은 시스템 제약(poppler-utils 미설치)으로 일부 논문에 대해 수행되지 못했으며, 상세한 정량 결과는 각 논문 원문을 참조하세요.*
