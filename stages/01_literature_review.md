# Stage 01 — Literature Review

## Workflow
1. `docs/01_literature_review.md` 가이드 확인
2. 문헌 검색 및 분석 수행
3. 산출물 → `outputs/literature_review/run_XX/`에 저장
4. 아래 지식 업데이트 (검증된 인사이트만 통합)
5. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 문헌 5편 이상 리뷰 완료
- SOTA 베이스라인 수치 최소 1개 확보
- 사용 가능한 데이터셋 확인 완료
- 기존 방법론의 한계(gap) 명시 완료

---

## 검증된 핵심 지식

### 섭동 예측 SOTA
| 모델 | R2 (sci-Plex) | 특징 |
|------|--------------|------|
| PerturbNet | 0.984 | 최고 성능 |
| PRnet | 0.969 | |
| CPA | 0.81 | 약물 임베딩에 MoA 클러스터링 관찰 (정량화 안 됨) |
| chemCPA (pretrained) | 0.68 (DEGs) | L1000 pretraining 2배 향상, 가장 직접적 비교 대상 |

### 대조학습 이론
- **Wang & Isola 2020**: alignment(positive pair 근접) + uniformity(초구면 균일 분포) = 대조 학습의 본질
  - L_align = E[||f(x)-f(x+)||^2], L_uniform = log E[exp(-t||f(x)-f(x')||^2)]
- **MoCL 2021**: 분자 그래프 대조 학습에 도메인 지식 주입 — bioisostere substitution(의미 보존 증강), global-level contrast

### MoA 분류 베이스라인
- **PANACEA DREAM Challenge**: 32 kinase inhibitor, 1,300 타겟, 21팀 — RNA-seq로 타겟 예측 검증
- **GPAR**: LINCS L1000, 978 gene → 103 MoA binary classifier, DNN > GSEA
- **DeepCE**: GNN+attention, L1000 1.4M 프로파일, de novo 화합물 발현 예측

### 대사 네트워크 최적화 관련
- **GEM→그래프 표현**: metabolite, reaction, gene 3노드타입 이종그래프 (HeteroData)가 자연적 표현
- **GNN 대리 모델**: HGTConv 2층 이종그래프 신경망이 대사 네트워크 임베딩에 적합 (PyTorch Geometric 지원)
- **dFBA 수치해석**: Euler 방법은 stiff system에서 불안정 (5.7x 과대추정 관측), BDF/Radau implicit solver 필수
- **FLYCOP**: 저장소 삭제(404), 사용 불가 → TOPSIS + Entropy weight로 객관화 대체 검증됨
- **COMETS**: v2.12.4 활성, Java 기반 → scipy BDF/Radau로 dFBA 직접 구현 가능

### 데이터
- **sci-Plex**: 188 compounds × 3 cell lines × 4 dosages = 649,340 cells, 7,561 genes, 16개 MoA 카테고리
- **COBRApy 모델**: textbook(95/72/137), iJO1366(2583/1805/1367), BiGG 108개 공개 모델

### Gap (다음 단계에 전달)
1. 섭동 예측 모델들은 DEGs R2만 최적화할 뿐, MoA 구조를 정량적으로 평가하지 않음 → MoA 클러스터링 품질을 명시적 타겟으로 하는 연구 부재
2. 대조 학습은 이미지/분자에서 검증되었으나, 단일세포 섭동 임베딩에 적용된 사례 없음
3. 대사 네트워크 최적화에서 FBA 기반 탐색은 조합 폭발에 직면하며, GNN 대리 모델 + Active Learning로 탐색 효율을 개선할 여지 존재하나 실증 부재
4. dFBA 기반 미생물 군집 시뮬레이션의 수치 안정성 문제가 상식적으로 알려져 있으나, 다목적 최적화(NSGA-II)와의 통합 사례 부재

---

## Run 이력 (세부 내용은 outputs/literature_review/run_XX/ 참조)
- run_01: LLM 바이오인포마틱스 10편 + Perturb-seq MoA 6모델. 갭 식별 성공, MoA 베이스라인/대조학습/데이터 미검증
- run_02: 대조학습 2편 + 섭동예측 3편 + MoA분류 3편 리뷰. 베이스라인/이론/데이터 기본 확보
