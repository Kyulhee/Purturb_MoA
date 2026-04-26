# Framing Document — Run 01

**Date**: 2026-04-26
**Stage**: 02 Framing
**Based on**: Literature Review run_01 + run_02

---

## 1. Research Question

> **대조 학습 기반 약물 임베딩이 sci-Plex 단일 세포 섭동 데이터에서 MoA(Mechanism of Action) 클러스터링 품질을 개선하는가?**

### 하위 질문
1. MoA-aware contrastive loss가 reconstruction-only latent space 대비 MoA 분리도를 얼마나 개선하는가?
2. 약물 분자 임베딩(rFCFP)과 발현 프로파일 임베딩의 정렬이 MoA 예측에 어떤 역할을 하는가?
3. leave-compound-out 설정에서 같은 MoA의 미지의 약물을 분류하는 정확도는?

---

## 2. Dataset

### Primary: sci-Plex3 (Srivatsan et al. 2020)
- **GEO Accession**: GSE139944 (sample: GSM4150378)
- **규모**: 188 compounds × 3 cell lines × 4 dosages = 649,340 cells, 7,561 genes
- **MoA 라벨**: `pathway_level_1` 기준 16개 카테고리 (원 논문 코드 확인)
  - Antioxidant, Apoptotic regulation, Cell cycle regulation, DNA damage & DNA repair, Epigenetic regulation, Focal adhesion signaling, HIF signaling, JAK/STAT signaling, Metabolic regulation, Neuronal signaling, Nuclear receptor signaling, PKC signaling, Protein folding & Protein degradation, TGF/BMP signaling, Tyrosine kinase signaling, Other
- **세포주**: A549 (폐암), MCF7 (유방암), K562 (백혈병)
- **용량**: 10nM, 100nM, 1µM, 10µM
- **데이터 포맷**: h5ad (chemCPA 전처리 버전 사용 가능)

### MoA 클래스 분포 (추정, 실제 검증 필요)
- 188 화합물 / 16 카테고리 → 평균 ~12 화합물/카테고리
- Epigenetic regulation이 최대 클래스 예상 (HDAC inhibitor 다수)
- "Other" 카테고리의 처리 방침 필요 (제외 vs 포함)
- **실제 분포는 데이터 다운로드 후 검증 필수**

### Pretraining: LINCS L1000 (선택적)
- ~1.4M 프로파일, ~20,000 화합물, ~50 세포주
- sci-Plex와 ~150개 화합물 + 모든 세포주 오버랩
- chemCPA에서 L1000 pretraining이 2배 성능 향상 검증

---

## 3. Evaluation Strategy

### 3.1 Primary Evaluation: Leave-Compound-Out Classification
- **설정**: 학습에 사용되지 않은 화합물의 MoA 분류
- **지표**: Top-1 Accuracy, F1 macro, AUROC (one-vs-rest)
- **의미**: "같은 MoA의 다른 화합물을 본 적이 있을 때, 새 화합물의 MoA를 맞출 수 있는가?"

### 3.2 Critical Evaluation: Leave-MoA-Out Clustering
- **설정**: 학습에 사용되지 않은 MoA 카테고리의 화합물이 latent space에서 클러스터링되는가
- **지표**: Silhouette score, ARI (Adjusted Rand Index), NMI (Normalized Mutual Information)
- **의미**: "완전히 새로운 MoA의 화합물들이 서로 유사한 임베딩을 갖는가?"
- **run_01 치명적 결함 교정**: 이전에는 leave-MoA-out을 분류 accuracy로 평가 → 구조적 불가능. 올바른 평가는 클러스터링 품질

### 3.3 Representation Quality Metrics
- **Alignment score**: 같은 MoA 약물 임베딩의 평균 거리 (Wang & Isola 기반)
  - `L_align = E[||f(x) - f(x+)||^2]` (같은 MoA pair)
- **Uniformity score**: 다른 MoA 약물 임베딩의 분포 균일성
  - `L_uniform = log E[exp(-t * ||f(x) - f(x')||^2)]` (다른 MoA pair)
- **MoA separation**: 클래스 간 중심 거리 vs 클래스 내 분산 비율

### 3.4 보조 평가: Gene Expression Reconstruction
- **지표**: r² (DEGs), MSE
- **의미**: MoA 분리가 재구성 품질을 희생하지 않는지 확인

---

## 4. Baselines (Quantitative)

### 계층 1: Random & Simple
| Baseline | 설정 | 예상 성능 | 근거 |
|----------|------|-----------|------|
| Random classifier | 16-class uniform | Top-1: 6.25% (1/16) | — |
| Majority class | 항상 최대 클래스 예측 | Top-1: ~15-20% 추정 | Epigenetic regulation 비율 추정 |
| GPAR-style DNN | 978 gene → MoA binary | AUROC > GSEA | GPAR 2021 (LINCS L1000) |

### 계층 2: Perturbation Prediction Models
| Baseline | 설정 | 성능 | 근거 |
|----------|------|------|------|
| CPA | 약물 임베딩만 사용 | MoA 클러스터링 관찰 (정량 미측정) | CPA 2023 |
| chemCPA (pretrained) | L1000→sci-Plex transfer | r²_DEGs = 0.68 | chemCPA 2022, sci-Plex3 leave-drug-covariate-out |

### 계층 3: MoA Classification Methods
| Baseline | 설정 | 성능 | 근거 |
|----------|------|------|------|
| PANACEA top methods | RNA-seq → target prediction | 21팀 성능 (AUROC 미추출) | PANACEA DREAM Challenge 2022 |
| DeepCE | GNN → expression → MoA | L1000 SOTA | DeepCE 2021 |

### 베이스라인 구현 우선순위
1. **Random + Majority class** (즉시 계산 가능)
2. **Simple DNN classifier** (sci-Plex 7,561 gene → 16-class) — GPAR 아키텍처 적용
3. **chemCPA drug embedding** (기존 코드 활용 가능 시)
4. **CMap connectivity** (참고용, 정량 수치 미확보)

---

## 5. Target Performance

### Leave-Compound-Out Classification
| 지표 | Random | Simple DNN (예상) | **Target** | 근거 |
|------|--------|-------------------|------------|------|
| Top-1 Accuracy | 6.25% | 40-60% | **>60%** | GPAR이 103 MoA에서 GSEA 대비 우수; 16-class에서는 더 높은 성능 기대 |
| F1 macro | ~0.06 | 0.35-0.50 | **>0.50** | 클래스 불균형 고려 |
| AUROC (OvR) | 0.50 | 0.85-0.90 | **>0.90** | 이진 분류 기준치 |

### Leave-MoA-Out Clustering (핵심 차별화 지표)
| 지표 | Reconstruction-only latent | **Target** | 근거 |
|------|---------------------------|------------|------|
| Silhouette score | ~0.05-0.10 (추정) | **>0.25** | CPA에서 MoA 클러스터링은 관찰되었으나 약함 |
| ARI | ~0.05 (random) | **>0.30** | 의미 있는 클러스터링의 최소 기준 |
| NMI | ~0.10 (random) | **>0.35** | 정보 이론적 독립성 기준 |

### Alignment & Uniformity
| 지표 | Reconstruction-only | **Target** | 근거 |
|------|---------------------|------------|------|
| Alignment (↓) | 높음 (MoA 무시) | **50% 감소** | Wang & Isola 이론 |
| Uniformity (적정) | 낮음 (재구성에 편향) | **개선** | 초구면 균일 분포 |

### 타겟 설정 원칙
- **보수적 타겟**: 베이스라인 대비 통계적으로 유의미(p < 0.05)한 개선
- **핵심 타겟**: leave-MoA-out ARI > 0.30 (이것이 달성되면 연구 목표 달성)
- **보조 타겟**: leave-compound-out accuracy > 60%, reconstruction r² > 0.50

---

## 6. Success Criteria

### Minimum Viable Result (MVR)
- leave-compound-out Top-1 accuracy > simple DNN baseline
- leave-MoA-out ARI > random baseline (유의미한 차이)
- Reconstruction r² > 0.40 (재구성 품질 유지)

### Target Result
- leave-compound-out Top-1 accuracy > 60%
- leave-MoA-out ARI > 0.30
- Alignment score 50% 개선 vs reconstruction-only

### Stretch Result
- leave-compound-out accuracy > PANACEA top methods
- leave-MoA-out ARI > 0.50
- L1000 pretraining으로 추가 성능 향상 검증

---

## 7. Data Verification Checklist

- [x] GEO accession 확인: GSE139944
- [x] sciPlex3 샘플 확인: GSM4150378
- [x] MoA 카테고리 16개 확인 (pathway_level_1)
- [ ] **실제 데이터 다운로드 및 로드** (진행 중)
- [ ] **MoA 라벨 분포 정량 확인** (클래스당 화합물 수)
- [ ] **클래스 불균형 분석** (min/max 클래스 비율)
- [ ] **유전자 발현 품질 확인** (missing rate, batch effect)
- [ ] **leave-compound-out / leave-MoA-out 분할 가능성 확인**

---

## 8. Risks & Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| 데이터 다운로드 실패 | Critical | GEO 직접 다운로드, chemCPA 전처리 버전 대안 |
| 16개 MoA 클래스 불균형 심각 | High | 소수 클래스 병합 또는 가중 손실 함수 |
| Leave-MoA-out에서 너무 적은 클래스 수 | Medium | 5-fold cross-validation, 클래스당 최소 3화합물 필터 |
| CMap 베이스라인 수치 미확보 | Low | Simple DNN 베이스라인으로 대체, CMap은 참고용 |
| Reconstruction-MoA 트레이드오프 | High | Loss weight 스케줄링, 단계적 학습 전략 |

---

## 9. Next Steps (to Planning)

1. **데이터 검증 완료**: MoA 분포, 클래스 불균형, 분할 가능성
2. **Simple DNN baseline 구현**: 7,561 gene → 16-class classifier로 첫 베이스라인 수치 확보
3. **Loss weight 설계**: L_recon vs L_contrastive 균형 (run_01에서 50:1 불균형이 치명적이었음)
4. **Encoder 동결 전략 재설계**: drug encoder fine-tuning 허용
5. **평가 파이프라인 구축**: leave-compound-out + leave-MoA-out + clustering metrics
