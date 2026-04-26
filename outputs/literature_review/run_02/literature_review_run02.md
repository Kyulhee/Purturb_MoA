# Literature Review Run 02 — Perturb-seq MoA Prediction

**Date**: 2026-04-26
**Scope**: run_01 결함 보완 — 대조학습, MoA 분류 베이스라인, CMap 검증, 데이터 검증

---

## 1. Contrastive Learning & Representation Quality

### 1.1 Wang & Isola (2020) — Alignment and Uniformity on the Hypersphere
- **Reference**: Wang & Isola, ICML 2020. arXiv:2005.10242
- **Key contribution**: 대조 표현 학습의 성공을 2가지 속성으로 분해:
  - **Alignment**: positive pair의 feature 유사도 (가까울수록 좋음)
  - **Uniformity**: feature 분포가 초구면 위에서 균일할수록 정보 보존이 좋음
- **핵심 공식**:
  - `L_align = E[||f(x) - f(x+)||^2]` (positive pair의 거리 최소화)
  - `L_uniform = log E[exp(-t * ||f(x) - f(x')||^2)]` (모든 sample pair의 거리 균일화, Gaussian kernel 기반)
  - `L_contrastive → L_align + λ * L_uniform` (M→∞ 극한에서)
- **정량 결과**: Alignment + Uniformity 직접 최적화가 SimCLR/MoCo와 동등 또는 우수한 downstream 성능 달성
- **프로젝트 적용**:
  - 약물 임베딩 공간에서 **같은 MoA의 약물은 alignment** (positive pair), **다른 MoA는 uniformity** (negative pair)로 학습
  - MoA-aware contrastive loss 설계의 이론적 근거 제공
  - L_recon과 L_contrastive의 gradient 균형 설계에 uniformity metric 활용 가능

### 1.2 MoCL (Sun et al. 2021) — Knowledge-aware Contrastive Learning from Molecular Graph
- **Reference**: Sun et al., KDD 2021. arXiv:2106.04509
- **Key contribution**: 분자 그래프 대조 학습에 도메인 지식 주입
  - **Local-level**: bioisostere substitution (의미 보존 증강) — 일반 증강(노드 삭제, 엣지 변경)은 분자 의미 파괴 가능
  - **Global-level**: 분자 간 유사도 정보를 contrastive scheme에 반영 (double contrast objective)
- **핵심 인사이트**:
  - 일반 그래프 증강(랜덤 노드/엣지 삭제)은 분자 그래프에 부적절 — 아스피린 phenyl ring의 탄소 원자 삭제는 방향성 시스템 파괴
  - "분자 특성을 바꾸지 않으면서 변화를 주는" 증강이 필요
- **정량 결과**: 다양한 분자 데이터셋에서 SOTA 달성 (linear + semi-supervised)
- **프로젝트 적용**:
  - rFCFP fingerprint 증강 전략 설계 시 참고 — 무작위 비트 flipping 대신 의미 보존 증강 필요
  - Global-level contrast: MoA 라벨 유사도를 대조 학습에 직접 반영하는 방식과 유사

---

## 2. Perturbation Prediction Models

### 2.1 chemCPA (Hetzel et al. 2022) — Predicting Cellular Responses to Novel Drug Perturbations
- **Reference**: Hetzel et al., NeurIPS 2022. arXiv:2204.13545
- **Key contribution**: CPA를 확장하여 **보이지 않는 약물(unseen drugs)** 에 대한 예측 가능
  - 분자 그래프 인코더(RDKit)를 통한 약물 임베딩 → perturbation latent space로 매핑
  - L1000 bulk 데이터로 pretraining → sci-Plex single-cell 데이터로 fine-tuning (transfer learning)
  - Architecture surgery: 유전자 수가 다른 데이터셋 간 전이 학습 지원
- **핵심 결과 (sci-Plex3 데이터, leave-drug-covariate-out)**:

| Dose | Model | E[r²]_all | E[r²]_DEGs | Median r²_all | Median r²_DEGs |
|------|-------|-----------|------------|---------------|----------------|
| 1µM | CPA | 0.72 | 0.54 | 0.86 | 0.67 |
| 1µM | chemCPA | 0.74 | 0.60 | 0.86 | 0.66 |
| 1µM | chemCPA_pretrained | **0.77** | **0.68** | 0.85 | **0.76** |
| 10µM | CPA | 0.54 | 0.34 | 0.52 | 0.26 |
| 10µM | chemCPA | 0.71 | 0.58 | 0.77 | 0.64 |
| 10µM | chemCPA_pretrained | **0.76** | **0.68** | **0.82** | **0.79** |

- **핵심 인사이트**:
  - L1000 pretraining이 성능 대폭 향상 (10µM DEGs: 0.34 → 0.68, 2배 개선)
  - Disentanglement: adversarial classifier로 약물/공변량 정보로부터 basal state 분리
  - 약물 인코더로 RDKit 분자 임베딩 사용 → **보이지 않는 약물에 대한 예측 가능**
- **프로젝트 적용**:
  - 우리 모델과 가장 직접적으로 비교 가능 — 같은 sci-Plex 데이터, 같은 leave-drug-out 평가
  - r²=0.68 (DEGs, pretrained)가 **발현 예측 베이스라인**
  - 단, 발현 예측 ≠ MoA 분류 — 우리 태스크는 클러스터링/분류 품질

### 2.2 CPA (Lotfollahi et al. 2023) — Compositional Perturbation Autoencoder
- **Reference**: Lotfollahi et al., Molecular Systems Biology, 2023. DOI:10.15252/msb.202211517. 228 citations
- **Key contribution**: 단일 세포 수준에서 섭동 응답 예측을 위한 해석 가능한 딥러닝
  - 선형 모델의 해석성 + 딥러닝의 유연성 결합
  - 보이지 않는 dosage, cell type, time point, species에 대한 일반화
  - 약물 조합 예측 검증 (baseline 대비 우수)
  - 5,329개 누락된 유전적 조합 imputation (97.6%)
- **프로젝트 적용**:
  - CPA는 약물 임베딩의 MoA 클러스터링이 관찰되었으나 정량화되지 않음 (run_01에서 이미 식별)
  - chemCPA가 CPA의 확장이므로, chemCPA 결과가 더 직접적 베이스라인

### 2.3 Jiang et al. (2022) — Therapeutic Algebra of Drug Responses at Single-Cell Resolution
- **Reference**: Jiang et al., 2022. arXiv:2208.10661
- **Key contribution**: 502개 면역조절 약물에 대한 단일 세포 반응 프로파일링 + 수학적 모델
  - 150만 개의 인간 면역 세포 프로파일링
  - 단일 규제 네트워크로 약물 단독 + 조합 반응 정량적 예측
  - **Drug response algebra**: 약물 조합의 가산적/비가산적 상호작용 분리
  - 약물 용량 titration을 통한 면역 세포 집단의 연속적 변조 예측
- **핵심 인사이트**:
  - 약물 조합 응답 = 선형 결합(linear) + 비선형 상호작용(non-additive)
  - 비가산적 응답은 개별 약물이 타겟하는 pathway 간 정보 통합 지점을 나타냄
- **프로젝트 적용**:
  - 약물 반응의 "대수(algebra)" 개념은 MoA 분류에 유용 — 같은 MoA 약물은 같은 유전자 발현 프로그램을 활성화해야 함
  - 단, 이 연구는 면역 세포 PBMC에 초점 — sci-Plex 암 세포와 다름

---

## 3. MoA Classification from Gene Expression — 베이스라인 성능

### 3.1 PANACEA DREAM Challenge (Douglass et al. 2022)
- **Reference**: Douglass et al., Cell Reports Medicine, 2022. DOI:10.1016/j.xcrm.2021.100492. 49 citations
- **Key contribution**: 최초의 대규모 MoA 예측 커뮤니티 챌린지
  - PANACEA: 25개 세포주 + ~400개 임상 항암제 + RNA-seq + dose-response
  - 32개 kinase inhibitor에 대해 1,300개 DrugBank 타겟 중 고친화도 결합 타겟 예측
  - 21개 팀 참가, compound identity blind
- **챌린지 설계**:
  - 입력: 약물-섭동된 RNA-seq 프로파일 + dose-response 곡선
  - 출력: 각 compound의 고친화도 결합 타겟 예측 (~1,300개 타겟 중)
  - 평가: proteome-wide polypharmacology 예측 정확도
- **우승 방법**: gene expression profile similarity analysis + deep learning (개별 데이터셋 학습)
- **핵심 인사이트**:
  - "Drug-perturbed RNA-seq data can be used to identify drug targets" (검증됨)
  - "Technology-based drug-target definitions often subsume literature definitions"
  - "Literature and screening datasets provide complementary information on drug mechanisms"
- **프로젝트 적용**:
  - **최초의 정량적 MoA 예측 베이스라인 제공** — 21팀의 성능이 비교 기준
  - 다만 kinase inhibitor에 국한 + target prediction (MoA classification과 다름)
  - 우리 태스크(MoA class 분류)와 가장 가까운 기존 벤치마크

### 3.2 DeepCE (Pham et al. 2021) — Mechanism-driven Phenotype Compound Screening
- **Reference**: Pham et al., Nature Machine Intelligence, 2021. DOI:10.1038/s42256-020-00285-9. 179 citations
- **Key contribution**: GNN + multihead attention으로 de novo 화합물의 발현 프로파일 예측
  - L1000 데이터의 불량 실험에서 유용한 정보 추출 (data augmentation)
  - 화학적 하위구조-유전자, 유전자-유전자 연관성 모델링
  - COVID-19 drug repurposing에 적용
- **데이터**: L1000 — ~1.4M 프로파일, ~50 세포주, ~20,000 화합물
- **CMap**: 원래 CMap은 ~1,300 화합물, 5개 암 세포주 (L1000은 확장판)
- **핵심 결과**: DeepCE > SOTA on L1000 gene expression prediction (de novo chemical setting)
  - 예측된 발현 프로파일의 downstream MoA 분류 효과도 검증됨
- **프로젝트 적용**:
  - L1000 → sci-Plex 전이 학습의 베이스라인 (chemCPA와 유사하지만 다른 접근)
  - 발현 예측 → MoA 분류의 2-stage 파이프라인 가능성 시사

### 3.3 GPAR (Gao et al. 2021) — Deep Learning for MOA from Gene Expression
- **Reference**: Gao et al., BMC Bioinformatics, 2021. DOI:10.1186/s12859-020-03915-6. 32 citations
- **Key contribution**: LINCS L1000 기반 MoA 예측을 위한 DNN 플랫폼
  - 978개 Landmark gene z-score를 feature로 사용
  - Binary classifier (positive: 같은 MoA 약물, negative: 낮은 transcript activity score + MoA 미주석)
  - 103개 MoA 평가 → **GPAR이 GSEA보다 우수**
  - DNN: 3 hidden layers (978→512→256), L1 정규화, RELU, dropout=0.1
- **핵심 인사이트**:
  - Deep learning의 hidden layer가 batch effect를 효과적으로 억제
  - DEG 기반 유사도(전통적 방법)는 batch effect와 common response에 취약
  - AUROC가 feature 수(유전자 수)에 비례하여 증가
- **프로젝트 적용**:
  - **LINCS 기반 MoA 분류의 직접적 베이스라인** — 978 gene → MoA binary classification
  - sci-Plex 7,561 gene → 더 많은 feature로 더 높은 성능 기대 가능
  - GPAR의 DNN 구조는 simple baseline으로 구현 가능

---

## 4. 데이터 검증 — sci-Plex

### 4.1 sci-Plex (Srivatsan et al. 2020)
- **Reference**: Srivatsan et al., Science, 2020. DOI:10.1126/science.aax6234. 355 citations
- **데이터 규모**:
  - 188 compounds × 3 cancer cell lines (A549, MCF7, K562) × 4 dosages (10nM, 100nM, 1µM, 10µM)
  - 649,340 cells, 7,561 drug-sensitive genes
  - **19개 MoA 클래스** (Srivatsan et al.가 모든 화합물에 할당)
- **데이터 접근**: GEO accession 확인 필요 (chemCPA 논문에서 사용 가능 확인)
- **L1000 오버랩**: ~150개 화합물 + 모든 세포주가 L1000과 겹침 (chemCPA에서 확인)
- **MoA 라벨**: 19개 MoA = 188개 화합물에 할당됨 (클래스당 평균 ~10개 화합물, 불균형 예상)
- **프로젝트 적용**:
  - **19개 MoA 클래스** — 분류 문제로서 적절한 규모
  - 클래스 불균형 분석 필요 — leave-MoA-out 시 일부 클래스의 샘플 수가 매우 적을 수 있음
  - L1000 오버랩 → pretraining 데이터 활용 가능 (chemCPA 전략)

---

## 5. 종합 분석 및 run_01 결함 해소 상태

| run_01 결함 | 해소 상태 | 근거 |
|------------|----------|------|
| MoA 분류 SOTA 부재 | **부분 해소** | PANACEA 챌린지(49 인용), GPAR(32 인용)가 베이스라인 제공. 그러나 정확한 수치(AUROC 등)는 PANACEA 논문에서 추가 추출 필요 |
| 대조학습 문헌 누락 | **해소** | Wang & Isola (alignment + uniformity), MoCL (분자 그래프 대조학습) 리뷰 완료 |
| 핵심 데이터셋 미검증 | **부분 해소** | sci-Plex: 188 화합물, 19 MoA, 649K cells 확인. GEO 다운로드는 실제 실행 필요 |
| 광범위 MoA 예측 문헌 누락 | **해소** | DeepCE, GPAR, PANACEA 리뷰 완료. LINCS → MoA 분류 파이프라인 확립 |
| CMap 정량적 성능 수치 | **미해소** | CMap은 전통적 bulk-level 기준이나, 정확한 MoA retrieval 수치(AUROC 등)를 제공하는 논문 미확보 |

---

## 6. 프로젝트를 위한 핵심 베이스라인 정리

### MoA 예측 베이스라인 계층
1. **Random baseline**: 19개 클래스 → 1/19 ≈ 5.3%
2. **GPAR-style DNN** (LINCS, 978 gene): MoA binary classifier — AUROC > GSEA
3. **PANACEA top methods** (RNA-seq, ~400 drugs): target prediction — 21팀 성능
4. **chemCPA** (sci-Plex, 7,561 gene): gene expression prediction r²=0.68 (DEGs), 약물 임베딩에 MoA 구조 관찰
5. **Our target**: MoA-aware contrastive learning으로 latent space에서 MoA 클러스터링 품질 개선

### 평가 지표 제안
- **Leave-compound-out classification**: Top-1 accuracy, F1 macro (기존과 동일)
- **Leave-MoA-out clustering**: Silhouette score, ARI, NMI (run_01에서 식별된 올바른 평가)
- **Alignment score**: 같은 MoA 약물 임베딩의 평균 거리 (Wang & Isola 기반)
- **Uniformity score**: 다른 MoA 약물 임베딩의 분포 균일성 (Wang & Isola 기반)

---

## 7. 다음 단계 (Framing으로 전달)

1. **데이터 검증**: sci-Plex 실제 다운로드 + MoA 라벨 분포(19개 클래스별 화합물 수) 확인
2. **베이스라인 수치화**: GPAR-style simple DNN을 sci-Plex에 적용하여 MoA 분류 AUROC 측정
3. **평가 전략 확정**: leave-MoA-out은 silhouette/ARI/NMI로 평가 (분류 accuracy 아님)
4. **연구 질문 명확화**: "대조 학습 기반 약물 임베딩이 sci-Plex 데이터에서 MoA 클러스터링 품질을 개선하는가?"
