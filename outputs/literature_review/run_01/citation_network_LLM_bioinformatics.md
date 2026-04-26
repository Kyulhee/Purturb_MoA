# 인용 네트워크 분석: LLM 기반 바이오인포매틱스 자동화

> 분석일: 2026-04-16 | 분석방법: OpenAlex CitationGraph (depth=2, both directions)

---

## 1. 인용 임팩트 랭킹

| 순위 | 논문 | 연도 | 인용수 | 바 |
|------|------|------|--------|----|
| 1 | GPTCelltype | 2024 | **205** | ████████████████████████████████████████ |
| 2 | AutoBA | 2024 | **51** | ██████████ |
| 3 | BixBench | 2025 | **49** | █████████ |
| 4 | DrBioRight 2.0 | 2025 | **33** | ██████ |
| 5 | BioResearcher | 2024 | **30** | ██████ |
| 6 | DrugAgent | 2024 | **36** | ███████ |
| 7 | Bioinfo-Bench | 2023 | **16** | ███ |
| 8 | PROTEUS | 2024 | **9** | █ |
| 9 | Agentic Bioinfo | 2025 | **6** | █ |
| 10 | ICGI | 2025 | **5** | █ |

> GPTCelltype가 압도적 허브 (205 인용). AutoBA와 BixBench가 2티어.

---

## 2. 연구 클러스터 분석

### Cluster A: 에이전트 기반 바이오인포 자동화
```
AutoBA ←→ Agentic Bioinfo ←→ BioResearcher
```
- **공통 특징**: 다중 에이전트 아키텍처, end-to-end 자동화, 도구 오케스트레이션
- **인용 교차**: AutoBA는 Agentic Bioinfo에 인용됨; BioResearcher도 동일 생태계
- **성장 동향**: 2025-2026년 ToolsGenie 2.0, BioAgents, STAnalyzer, DynaMate 등 급속 확장

### Cluster B: Single-Cell / 오믹스 LLM 응용
```
GPTCelltype ←→ DrBioRight 2.0 ←→ PROTEUS
```
- **공통 특징**: 도메인 특화 LLM 사용, 데이터 기반 주석 달기, 프로테오믹스/유전체학
- **인용 교차**: GPTCelltype이 205 인용으로 이 클러스터의 핵심 허브
- **성장 동향**: 성숙기 — GPTCelltype을 기반으로 ReCellTy, CellEmbed 등 후속 연구 활발

### Cluster C: LLM 평가 & 벤치마킹
```
Bioinfo-Bench ←→ BixBench
```
- **공통 특징**: LLM 능력 평가, 한계 노출, 프레임워크 설계
- **인용 교차**: 두 벤치마크 모두 "LLM이 실무에서 아직 부족하다"는 결론으로 수렴
- **성장 동향**: 중요 갭 — BixBench에서 프론티어 모델 17% 정답률. 2026년 Genome Biology에 single-cell omics 벤치마크 게재

### Cluster D: 신약 발견 & 인과 추론
```
DrugAgent ←→ ICGI
```
- **공통 특징**: 제약 ML 다중 에이전트, 인과 프롬프팅, 유전자 식별
- **인용 교차**: "Next-generation agentic AI for healthcare" (72 인용)가 두 논문 모두 인용
- **성장 동향**: 신흥 분야 — 2025-2026년 차세대 에이전틱 AI가 헬스케어로 확장

---

## 3. 브릿지 논문 (여러 클러스터를 연결)

| 브릿지 논문 | 연도 | 인용 | 연결하는 논문 |
|-------------|------|------|--------------|
| Bioinformatics with ChatGPT: Year One Review | 2024 | 29 | AutoBA, GPTCelltype |
| The potential of LLMs to advance precision oncology | 2025 | 18 | DrBioRight 2.0, ICGI |
| AI Agents vs. Agentic AI: Conceptual Taxonomy | 2025 | 62 | DrugAgent |
| Next-generation agentic AI for healthcare | 2025 | 72 | DrugAgent, ICGI |
| Prompt-based bioinformatics: new interface | 2025 | 3 | DrBioRight 2.0 |
| General-Purpose Models for Chemical Sciences | 2026 | 2 | DrugAgent, BixBench |
| AI agents for biological research: a survey | 2026 | 1 | Bioinfo-Bench |

---

## 4. 2025-2026 신흥 트렌드 (인용 네트워크에서 발견)

### 4.1 에이전틱 AI의 헬스케어 확장
- "Next-generation agentic AI for transforming healthcare" (72 인용) — DrugAgent 인용
- "AI Agents vs. Agentic AI: A Conceptual Taxonomy" (62 인용) — 에이전트 vs 에이전틱 AI 구분
- 방향: 단순 도구 호출 → 자율적 의사결정 시스템으로 진화

### 4.2 도메인 특화 다중 에이전트 시스템 급증
- **ToolsGenie 2.0** (2026) — 확장 가능한 바이오인포 자동화 멀티에이전트 (BixBench 인용)
- **BioAgents** (2025, Nat. Sci. Rep.) — 바이오인포 분석 멀티에이전트 (Bioinfo-Bench 인용)
- **STAnalyzer** (2026) — 공간 전사체학 에이전틱 아키텍처
- **ChatMolData** (2025) — 분자 데이터 처리 멀티모달 에이전트
- **DynaMate** (2025) — 맞춤형 연구 워크플로우 AI 에이전트

### 4.3 LLM 에이전트의 Single-Cell 오믹스 벤치마킹
- "Benchmarking LLM-based agents for single-cell omics analysis" (Genome Biology 2026, 2 인용)
- BixBench의 한계를 single-cell 도메인에서 심화 평가

### 4.4 자가 진화형 자율 연구 시스템
- "Empowering AI data scientists using multi-agent LLM framework with self-evolving capabilities" (Nature BME 2026, 0 인용)
- BioResearcher의 end-to-end 자동화를 넘어 **자가 학습/진화** 단계로

### 4.5 인과 추론 + LLM 결합
- ICGI의 인과 프롬프팅 접근이 암 유전자 식별을 넘어 다른 오믹스 도메인으로 확산
- "Decoding critical targets using LLMs" (2024, 4 인용) — EBV 질환에 유사 적용

### 4.6 프롬프트 기반 바이오인포매틱스 — 새 패러다임
- "Prompt-based bioinformatics: a new interface for multi-omics analysis" (Nature Reviews Genetics 2025, 3 인용)
- DrBioRight 2.0 인용 — 자연어가 바이오인포의 새 인터페이스가 됨

### 4.7 인간-AI 협업 프레임워크
- "A conceptual framework for human–AI collaborative genome annotation" (Brief. Bioinf. 2025)
- 완전 자율이 아닌 **인간-인-더-루프** 패러다임 정착 시그널

---

## 5. 핵심 인사이트

1. **GPTCelltype (205 인용)이 압도적 허브** — 특정 태스크(세포타입 주석)에 LLM을 정밀 적용하는 전략이 가장 높은 임팩트. 범용 자동화보다 도메인 특화가 인용에서 유리.

2. **AutoBA와 Agentic Bioinformatics가 에이전트 자동화 클러스터의 핵심** — 2025-2026년 후속 시스템(ToolsGenie 2.0, BioAgents 등)이 이 논문들을 기반으로 분화.

3. **BixBench(17% 정답률)가 현실 점검 역할** — 낮은 성능 수치가 오히려 다음 세대 시스템 개발의 자극제. 후속 벤치마크들이 BixBench를 기준선으로 삼음.

4. **필드가 3방향으로 분화 중**:
   - (a) **도메인 특화 도구** (GPTCelltype, DrBioRight, PROTEUS)
   - (b) **범용 에이전트 프레임워크** (AutoBA, BioResearcher, DrugAgent)
   - (c) **평가/벤치마킹** (BixBench, Bioinfo-Bench)

5. **2025-2026 수렴 트렌드**: 에이전트 프레임워크 + 도메인 지식 + 자가 진화 능력의 결합. "Self-evolving multi-agent LLM" (Nature BME 2026)이 다음 단계의 청사진.

6. **인간-AI 협업이 현실적 방향** — 완전 자율(BixBench 17%)은 요원; 인간-인-더-루프 프레임워크가 실용적 합의점.

7. **가장 주목할 신흥 논문**:
   - "Prompt-based bioinformatics" (Nat. Rev. Genet. 2025) — 새 패러다임 제시
   - "Pushing the boundaries of autonomous biological discovery" (Nat. Methods 2026) — Agentic Bioinfo 인용
   - "Empowering AI data scientists with self-evolving capabilities" (Nat. BME 2026) — BixBench 인용

---

## 6. 네트워크 토폴로지 요약

```
                    ┌──────────────────┐
          ┌────────►│  GPTCelltype     │◄──── Single-Cell/오믹스 클러스터
          │    205  │  (205 cit, 허브)  │        (DrBioRight, PROTEUS)
          │         └──────────────────┘
          │
   ┌──────┴──────┐     ┌──────────────────┐
   │ Bioinfo-Bench│────►│    BixBench      │◄─── 벤치마킹 클러스터
   │ (16 cit)     │     │   (49 cit)       │
   └──────────────┘     └────────┬─────────┘
                                 │
          │                      │ 인용
          │              ┌───────▼──────────┐
   ┌──────┴──────┐       │  AutoBA          │◄─── 에이전트 자동화 클러스터
   │ DrugAgent   │──────►│  (51 cit)        │     (BioResearcher, Agentic Bioinfo)
   │ (36 cit)    │       └──────────────────┘
   └──────┬──────┘
          │         ┌──────────────────┐
          └────────►│     ICGI         │◄─── 인과추론 클러스터
                    │   (5 cit)        │
                    └──────────────────┘
```

---

*이 분석은 OpenAlex CitationGraph API를 통해 수집된 2-hop 인용 네트워크 데이터를 기반으로 작성되었습니다.*
