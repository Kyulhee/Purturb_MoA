# 연구 플랜: Perturb-seq 데이터 기반 약물 MoA(Mechanism of Action) 예측

> 작성일: 2026-04-16 | 최종 업데이트: 2026-04-16 | 근거 수준: 6편 전문 파싱 + 3편 초록/메타데이터

---

## 0. 근거 수준 명시

| 항목 | 근거 수준 | 비고 |
|------|----------|------|
| 핵심 논문 6편 (GEARS, CPA, PerturbNet, PRnet, UNAGI, CellOT) | **전문 파싱 완료** | 아키텍처, loss, 벤치마크 수치 확인 |
| 데이터셋 논문 3편 (Replogle, Norman, sci-Plex) | **초록 수준** | 전문 다운로드 실패 (페이월); 타 논문 인용 수치 활용 |
| 제안 아키텍처 | **연구자 제안** (미검증) | 기존 모델 구성요소 기반이나, 통합 방식은 새로운 제안 |
| 공개 데이터셋 접근성 | **도구 검증 필요** | GEO 실제 다운로드 테스트 미수행 |
| MoA-aware contrastive loss | **연구자 제안** (미검증) | SimCLR/MoCo 문헌에서 검증된 loss family이나, cross-perturbation 적용은 새로움 |

---

## 1. 연구 배경 & 핵심 질문

### 1.1 배경
- Perturb-seq(CRISPR + scRNA-seq)는 유전자 perturbation이 세포 전사체에 미치는 영향을 대규모로 매핑
- 약물의 작용 기전(MoA)은 유전자 perturbation 효과와 유사한 전사체 변화를 유발하는 경우가 많음
- **핵심 아이디어**: CRISPR perturbation profile ↔ drug-induced signature 간 매칭 → 약물 MoA 예측
- 기존 Connectivity Map(CMap)은 bulk 수준이지만, single-cell 해상도에서의 MoA 예측은 아직 초기 단계

### 1.2 핵심 연구 질문
> **Perturb-seq 데이터에서 추출한 유전자 perturbation signature를 활용하여, 신약 후보 물질의 작용 기전(MoA)을 single-cell 해상도에서 얼마나 정확하게 예측할 수 있는가?**

### 1.3 세부 질문
1. CRISPR KO perturbation profile과 약물 처리 profile 간의 전사체 수준 유사도가 MoA 예측에 유효한 지표인가?
2. 기존 perturbation 예측 모델(GEARS, CPA, PerturbNet)이 학습한 latent space가 drug MoA 분류에 전이 학습 가능한가?
3. Single-cell 해상도의 perturbation signature가 bulk CMap 대비 MoA 예측에서 어떤 추가 정보를 제공하는가?

---

## 2. 선행연구 심층 분석 결과 (전문 파싱 기반)

### 2.1 확인된 연구 갭

> **핵심 갭**: CRISPR perturbation profile과 drug perturbation profile을 single-cell 해상도에서 직접 매칭하여 MoA를 예측하는 연구가 **존재하지 않음**

| 연구 방향 | 기존 연구 상태 | 갭 |
|-----------|---------------|-----|
| CRISPR perturbation 효과 예측 | GEARS, PerturbNet (genetic) | — |
| Drug perturbation 효과 예측 | CPA, PerturbNet (chemical), PRnet | — |
| **CRISPR ↔ Drug cross-perturbation** | **❌ 없음** | **이 연구의 타겟** |
| Drug MoA 분류/예측 | PRnet (GSEA), UNAGI (CMAP) | Bulk 수준에 머물름 |
| **Single-cell 해상도 MoA** | **❌ 없음** | **이 연구의 타겟** |
| Latent space MoA 클러스터링 | CPA (발견적), CellOT (계층적) | 분류 task로 발전시킨 연구 없음 |

### 2.2 핵심 모델 비교

| 모델 | 기본 구조 | 약물 처리 | MoA 직접 예측 | sci-Plex R² | 핵심 한계 |
|------|----------|----------|-------------|-------------|----------|
| **GEARS** | GNN + gene regulatory network | ❌ | ❌ | — | Drug 미지원; poorly connected gene 예측 불가 |
| **CPA** | VAE + adversarial disentanglement | ✅ (chemCPA) | ❌ (단, MoA 클러스터링 발견) | 0.81 | 클러스터링 → 분류 발전 없음 |
| **PerturbNet** | cINN + ChemicalVAE + GenotypeVAE | ✅ | ❌ | **0.984** | CRISPR-Drug cross-matching 없음 |
| **PRnet** | Deep generative + rFCFP 임베딩 | ✅ | ✅ (GSEA) | 0.969 | GSEA bulk 기반; CRISPR 연결 없음 |
| **UNAGI** | VAE-GAN + GCN | ✅ (CMAP) | ✅ (disease→drug) | — | CMAP bulk 한계; Perturb-seq 미활용 |
| **CellOT** | Neural OT (ICNN) | ✅ | ❌ (단, transport cost 군집화) | — | Per-perturbation 모델; 확장성 제한 |

### 2.3 데이터셋 명세

| 데이터셋 | Perturbation 유형 | 규모 | Cell line | MoA 관련성 | 접근성 |
|----------|-------------------|------|-----------|------------|--------|
| **Replogle 2022** | CRISPRi (genome-scale) | >2.5M cells, 5,000+ genes | K562 | ⭐⭐⭐ | GEO (추정) |
| **Norman 2019** | CRISPRa/CRISPRi (2-gene 조합) | ~218K cells, 131 조합 | K562 | ⭐⭐ | GEO 확인 필요 |
| **sci-Plex** | 약물 (188 compounds, 4 doses) | ~650K cells | 3 cancer cell lines | ⭐⭐⭐⭐⭐ | GEO (PerturbNet 전처리 코드 활용 가능) |
| **LINCS L1000** | 약물 (수만 compounds) | ~100M bulk | 88 cell lines | ⭐⭐⭐⭐ | 공개 (clue.io) |

> **근거 수준**: Replogle, Norman, sci-Plex은 초록 수준. sci-Plex 상세 수치는 PerturbNet/PRnet 전문에서 인용. LINCS L1000은 PRnet 전문에서 확인.

---

## 3. 연구 방법론: Cross-Perturbation MoA Framework

### 3.1 제안 아키텍처 개요

기존 모델 분석에서 도출한 핵심 설계 원칙:
- **PerturbNet으로부터**: cINN 양방향 매핑 + GenotypeVAE (GO-based) + ChemicalVAE
- **CPA로부터**: Perturbation effect의 가산적 분해 (disentanglement) + adversarial training
- **PRnet으로부터**: rFCFP 임베딩 (fingerprint × dose) + GSEA 검증 프레임워크
- **CellOT로부터**: 분포 수준 매핑의 필요성 (mean shift 불충분)
- **새로운 제안**: MoA-aware contrastive loss로 cross-perturbation latent space 학습

```
┌─────────────────────────────────────────────────────┐
│           Cross-Perturbation MoA Framework           │
├─────────────────────────────────────────────────────┤
│                                                      │
│  ┌──────────────┐         ┌──────────────┐           │
│  │  CRISPR KO   │         │    Drug      │           │
│  │  Perturbation│         │  Perturbation│           │
│  └──────┬───────┘         └──────┬───────┘           │
│         │                        │                    │
│  ┌──────▼───────┐         ┌──────▼───────┐           │
│  │ GenotypeVAE  │         │ Drug Encoder │           │
│  │ (GO-based)   │         │ (rFCFP+dose) │           │
│  └──────┬───────┘         └──────┬───────┘           │
│         │                        │                    │
│         └────────┬───────────────┘                    │
│                  │                                    │
│         ┌────────▼────────┐                           │
│         │  Shared Latent  │                           │
│         │     Space       │                           │
│         │  (Disentangled  │                           │
│         │   + MoA-aware)  │                           │
│         └────────┬────────┘                           │
│                  │                                    │
│     ┌────────────┼────────────┐                      │
│     │            │            │                      │
│  ┌──▼──┐   ┌────▼────┐  ┌───▼───┐                   │
│  │MoA  │   │Cell     │  │Cross  │                   │
│  │Class│   │State    │  │Perturb│                   │
│  │Head │   │Decoder  │  │Match  │                   │
│  └─────┘   └─────────┘  └───────┘                   │
│                                                      │
│  Loss: Recon + Adversarial + MoA-Contrastive         │
└─────────────────────────────────────────────────────┘
```

### 3.2 핵심 구성요소 상세

#### 3.2.1 GenotypeVAE (CRISPR perturbation encoder)

**근거**: PerturbNet의 GenotypeVAE가 GO annotation 기반으로 ~177M gene combination 학습 (전문 파싱 근거)

```
입력: CRISPRi/KO 대상 유전자 목록
인코딩:
  - Gene → GO annotation embedding (15,988 terms, 18,832 genes)
  - Perturbation effect = Σ(gene_i embedding) [가산적, CPA-style]
  - Basal state와 adversarial 분리
출력: z_perturb (perturbation latent vector)
```

**주의**: Replogle 2022가 CRISPRi (partial knockdown)를 사용 → KO보다 drug inhibition에 더 근접. 이는 "CRISPR KO ≠ drug inhibition" 리스크의 완화 요인.

#### 3.2.2 Drug Encoder (rFCFP-based)

**근거**: PRnet의 rFCFP 임베딩이 PCC=0.8 (unseen compounds) 달성 (전문 파싱 근거)

```
입력: SMILES 문자열 + 용량(dose)
인코딩:
  - SMILES → FCFP fingerprint (RDKit)
  - rFCFP = FCFP × log10(dose + 1)  [PRnet 방식]
  - MLP → z_drug (drug latent vector)
출력: z_drug
```

**선택 근거**: ChemicalVAE (PerturbNet) 대비 rFCFP 채택 이유:
1. 해석 가능성: fingerprint → 특정 chemical substructure 추적 가능
2. 용량 정보 내장: dose-response 곡선의 핵심 정보
3. 경량성: VAE 사전학습 (ZINC 250K) 불필요

#### 3.2.3 Shared Latent Space + MoA-aware Contrastive Loss

**근거**: CPA의 drug embedding이 MoA별 클러스터링됨 (전문 파싱 근거). 단, 이를 분류 task로 활용한 연구 없음 → 본 연구의 핵심 기여.

```
Contrastive Loss:
  - 같은 MoA 그룹: 약물 ↔ CRISPR perturbation을 latent space에서 가깝게
    L_attract = -log(exp(sim(z_i, z_j)/τ) / Σ exp(sim(z_i, z_k)/τ))
  
  - 다른 MoA 그룹: 멀게
    Hard negative: 같은 pathway 내 다른 MoA (예: HDAC1 KO vs CDK4 inhibitor)
  
  - Alignment + Uniformity (Wang & Isola 2022):
    L_align = ||z_i - z_j||²  (same MoA)
    L_uniform = log Σ exp(-2t||z_i - z_j||²)  (전체 분포 균일성)
```

**[미검증 — 연구자 제안]**: CRISPR KO와 drug perturbation이 같은 MoA에서 실제로 latent space에서 가까워지는지는 실험적으로 검증해야 함. CPA의 클러스터링 결과가 간접적 지지이나, cross-perturbation (CRISPR↔Drug) 설정에서의 검증은 없음.

#### 3.2.4 세 가지 출력 헤드

1. **MoA Classification Head**: 약물 perturbation의 MoA 라벨 분류 (ATC code 기반)
   - 학습: 약물 데이터 (sci-Plex 188 compounds)의 MoA 라벨로 supervised 학습
   - 추론: CRISPR KO의 perturbation profile → MoA 분류 (cross-perturbation)

2. **Cell State Decoder**: Perturbation 효과 재구성
   - 평가: R², Pearson correlation (예측 vs 실제)
   - 기능: Latent space 품질 보장 (재구성 가능해야 의미 있는 space)

3. **Cross-Perturbation Matching**: CRISPR KO ↔ Drug profile 유사도
   - 평가: Latent space에서 CRISPR perturbation과 가장 유사한 drug perturbation 검색
   - 기능: MoA 추론의 해석 가능성 제공 ("이 CRISPR KO는 HDAC inhibitor와 유사한 효과")

### 3.3 손실 함수

```
L_total = λ₁ · L_recon          # Cell state 재구성 (MSE/R²)
        + λ₂ · L_adversarial    # Perturbation ↔ covariate 분리 (CPA-style)
        + λ₃ · L_contrastive    # MoA-aware contrastive (본 연구 핵심)
        + λ₄ · L_classification # MoA 분류 (cross-entropy)
```

λ 하이퍼파라미터는 ablation study로 결정.

---

## 4. 평가 프로토콜

### 4.1 데이터 분할 전략

| 전략 | 설명 | 평가 목적 |
|------|------|----------|
| **Leave-compound-out** | 특정 약물을 test에 배제 | Unseen compound MoA 예측 |
| **Leave-MoA-out** | 특정 MoA 클래스 전체를 test에 배제 | Zero-shot MoA 분류 |
| **Leave-cell-line-out** | 특정 cell line을 test에 배제 | Cross-cell-line 일반화 |
| **Cross-perturbation** | CRISPR KO로 학습 → 약물 MoA 예측 | **핵심: CRISPR→Drug 전이** |

### 4.2 평가 지표

| 지표 | 용도 | Baseline 기준 |
|------|------|---------------|
| MoA 분류 정확도 (Top-1/3/5) | MoA 예측 성능 | Random, CMap, PRnet |
| R² (gene-level) | Perturbation 효과 재구성 | PerturbNet 0.984, PRnet 0.969 |
| Pearson correlation | Latent space 유사도 | PRnet PCC=0.8 |
| ARI / NMI | MoA 클러스터링 품질 | CPA 클러스터링 |
| AUROC (MoA binary) | 특정 MoA 탐지 | CMap connectivity |

### 4.3 Baseline 비교

| Baseline | 비교 포인트 |
|----------|------------|
| **CMap (bulk)** | Bulk 대비 single-cell 이점 |
| **CPA (원본)** | Drug embedding 클러스터링 |
| **PerturbNet** | Perturbation 예측 정확도 |
| **PRnet (GSEA)** | Drug screening 접근 |
| **Random Forest** | 단순 분류기 upper bound |

---

## 5. 실험 계획

### Phase 1: 데이터 준비 & Baseline 재현 (1-2주)

1. **sci-Plex 데이터 전처리**
   - PerturbNet/CPA/PRnet 공개 코드에서 전처리 파이프라인 재현
   - 648,737 cells, 5,087 genes, 180 drug treatments (PerturbNet 전처리 기준)
   - MoA annotation: DrugBank/ATC code 매핑

2. **Replogle 2022 (CRISPRi) 데이터 접근**
   - GEO accession 확인 후 다운로드
   - CRISPRi 데이터 = drug inhibition에 더 근접 → 리스크 완화

3. **CPA drug embedding MoA 클러스터링 재현**
   - 공식 코드로 sci-Plex 학습
   - Drug embedding t-SNE/UMAP + MoA별 ARI/NMI 정량화
   - **근거**: CPA 논문이 클러스터링 보고했으나 정량화 없음 → 본 연구의 baseline

### Phase 2: 모델 프로토타입 (2-3주)

1. **GenotypeVAE + Drug Encoder → Shared Latent Space**
   - PyTorch + PyG 구현
   - GenotypeVAE: GO annotation 기반 (PerturbNet 방식)
   - Drug Encoder: rFCFP (PRnet 방식)
   - CPA-style disentanglement + adversarial training

2. **MoA-aware Contrastive Learning 실험**
   - Temperature scaling, hard negative mining
   - Alignment + Uniformity loss (Wang & Isola 2022)
   - λ 하이퍼파라미터 탐색

3. **초기 평가**
   - Leave-compound-out MoA 분류 정확도
   - Latent space t-SNE/UMAP 시각화
   - Cross-perturbation matching 정성 평가

### Phase 3: 전체 실험 & 분석 (2-3주)

1. **전체 Baseline 비교**
2. **Ablation Study**: Contrastive loss 유무, latent space 차원, transfer learning 전략
3. **케이스 스터디**: HDAC inhibitor, CDK inhibitor 예측 성공/실패 분석
4. **한계 분석**: Cell line specificity, CRISPRi↔drug gap, off-target
5. **논문 초안**

---

## 6. 타임라인 요약

| 주차 | Phase | 산출물 |
|------|-------|--------|
| 1-2주 | Phase 1: 데이터 준비 & Baseline | 전처리 파이프라인, CPA 클러스터링 재현 결과 |
| 3-5주 | Phase 2: 모델 프로토타입 | Cross-Perturbation MoA Framework 프로토타입, 초기 결과 |
| 6-8주 | Phase 3: 전체 실험 & 분석 | 전체 평가, ablation, 논문 초안 |

---

## 7. 핵심 리스크 & 완화 전략 (전문 분석 기반 업데이트)

| 리스크 | 확률 | 심각도 | 완화 전략 | 근거 |
|--------|------|--------|----------|------|
| **CRISPR KO ≠ Drug inhibition** | 높음 | 높음 | Replogle 2022 CRISPRi 데이터 활용 (partial knockdown → drug에 더 근접); stoichiometry 보정 항 추가 | Replogle 초록에서 CRISPRi 확인 |
| **Latent space 정렬 불가** | 중간 | 높음 | Warm-start: 각 VAE 사전학습 후 공통 space로 미세조정 | PerturbNet은 별도 VAE → 공통 space 미검증 |
| **sci-Plex 3 cell line 일반화 부족** | 중간 | 중간 | LINCS L1000 (88 cell lines, bulk) 보조 학습 | PRnet이 88 cell lines 활용 (전문 근거) |
| **MoA 라벨 불균형** | 중간 | 중간 | Focal loss; few-shot 대비 metric learning | sci-Plex 188 compounds MoA 분포 확인 필요 |
| **Contrastive loss trivial solution** | 낮음–중간 | 중간 | Temperature scaling; hard negative mining; alignment+uniformity loss | SimCLR/MoCo에서 검증된 완화법 |
| **Novelty 확보** | 중간 | 높음 | Cross-perturbation contrastive learning = 핵심 차별화; CRISPR→Drug zero-shot MoA 선행 없음 | 본 분석에서 확인된 갭 |

---

## 8. 예상 논문 구조

1. **Introduction**: Perturb-seq + drug MoA 예측의 필요성, 기존 한계
2. **Related Work**: GEARS, CPA, PerturbNet, PRnet 분석 및 갭 식별
3. **Methods**: Cross-Perturbation MoA Framework
   - Shared latent space 설계
   - MoA-aware contrastive learning
   - Cross-perturbation inference
4. **Results**:
   - Baseline 대비 MoA 분류 성능
   - Zero-shot MoA 예측 (Leave-MoA-out)
   - CRISPR→Drug cross-perturbation 전이
   - Latent space 해석 (MoA 클러스터링, attention/saliency)
   - 케이스 스터디 (HDAC inhibitor, CDK inhibitor)
5. **Discussion**: CRISPR-drug 갭, 한계, 향후 방향
6. **Data & Code Availability**: 재현 가능한 파이프라인 공개

---

## 9. 즉각적 다음 단계

1. **sci-Plex 데이터 전처리 파이프라인 구축** — PerturbNet/CPA/PRnet 공개 코드 활용
2. **MoA annotation 매핑** — DrugBank/ATC code로 188 compounds 분류
3. **CPA baseline 재현** — Drug embedding MoA 클러스터링 정량화
4. **Replogle 2022 GEO accession 확인** — CRISPRi 데이터 접근

---

*이 플랜은 6편의 전문 PDF 파싱 결과와 3편의 초록/메타데이터를 기반으로 작성되었습니다. 제안 아키텍처는 기존 모델의 검증된 구성요소(GenotypeVAE, rFCFP, CPA disentanglement)를 통합하되, MoA-aware contrastive learning과 cross-perturbation matching은 [미검증 — 연구자 제안]입니다.*
