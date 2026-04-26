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

---

## 검증된 핵심 지식

### 섭동 예측 SOTA
| 모델 | R² (sci-Plex) | 특징 |
|------|--------------|------|
| PerturbNet | 0.984 | 최고 성능 |
| PRnet | 0.969 | |
| CPA | 0.81 | 약물 임베딩에 MoA 클러스터링 관찰 (정량화 안 됨) |
| chemCPA (pretrained) | 0.68 (DEGs) | L1000 pretraining 2배 향상, 가장 직접적 비교 대상 |

### 대조학습 이론
- **Wang & Isola 2020**: alignment(positive pair 근접) + uniformity(초구면 균일 분포) = 대조 학습의 본질
  - L_align = E[||f(x)-f(x+)||²], L_uniform = log E[exp(-t||f(x)-f(x')||²)]
- **MoCL 2021**: 분자 그래� 대조 학습에 도메인 지식 주입 — bioisostere substitution(의미 보존 증강), global-level contrast

### MoA 분류 베이스라인
- **PANACEA DREAM Challenge**: 32 kinase inhibitor, 1,300 타겟, 21팀 — RNA-seq로 타겟 예측 검증
- **GPAR**: LINCS L1000, 978 gene → 103 MoA binary classifier, DNN > GSEA
- **DeepCE**: GNN+attention, L1000 1.4M 프로파일, de novo 화합물 발현 예측

### 데이터: sci-Plex
- 188 compounds × 3 cell lines × 4 dosages = 649,340 cells, 7,561 genes
- **16개 MoA 카테고리** (pathway_level_1 기준, 원논문 코드에서 확인)
- L1000과 ~150개 화합물 오버랩 → pretraining 가능

### 미해결
- CMap MoA retrieval 정량 수치(AUROC 등) 미확보
- 실제 데이터 다운로드 및 MoA 분포 정량 확인 필요

### 다음 단계에 전달
1. sci-Plex 실제 다운로드 + MoA 라벨 분포 확인
2. leave-MoA-out은 silhouette/ARI/NMI로 평가 (분류 accuracy 아님)
3. 베이스라인 계층: random(6.25%) → GPAR-style DNN → PANACEA → chemCPA → our method

---

## Run 이력 (세부 내용은 outputs/literature_review/run_XX/ 참조)
- run_01: LLM 바이오인포마틱스 10편 + Perturb-seq MoA 6모델. 갭 식별 성공, MoA 베이스라인/대조학습/데이터 미검증
- run_02: 대조학습 2편 + 섭동예측 3편 + MoA분류 3편 리뷰. 베이스라인/이론/데이터 기본 확보
