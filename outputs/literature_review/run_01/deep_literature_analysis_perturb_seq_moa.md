# 심층 문헌 분석 보고서: Perturb-seq 기반 약물 MoA 예측

> 작성일: 2026-04-16 | 근거 수준: 6편 전문 파싱 + 3편 초록/메타데이터

---

## 0. 근거 수준 명시

| 항목 | 근거 수준 | 비고 |
|------|----------|------|
| GEARS | **전문 파싱 완료** | 아키텍처, loss, 평가 프로토콜 확인 |
| CPA | **전문 파싱 완료** | compositional decomposition, chemCPA 확장 확인 |
| PerturbNet | **전문 파싱 완료** | cINN, ChemicalVAE, GenotypeVAE, 벤치마크 수치 확인 |
| PRnet | **전문 파싱 완료** | rFCFP, GSEA 스크리닝, 실험 검증 수치 확인 |
| UNAGI | **전문 파싱 완료** | VAE-GAN+GCN, IPF 적용, CMAP 스크리닝 확인 |
| CellOT | **전문 파싱 완료** | Neural OT, ICNN, 일반화 실험 확인 |
| Replogle 2022 | **초록 수준** | 전문 다운로드 실패 (Cell 저널 페이월) |
| Norman 2019 | **초록 수준** | 전문 다운로드 실패 (Science 저널 페이월) |
| sci-Plex | **초록 수준** | 전문 다운로드 실패 (Science 저널 페이월), 타 논문에서 인용된 상세 수치 활용 가능 |

---

## 1. 핵심 모델 비교 분석

### 1.1 아키텍처 비교

| 차원 | GEARS | CPA | PerturbNet | PRnet | UNAGI | CellOT |
|------|-------|-----|------------|-------|-------|--------|
| **기본 구조** | GNN + cross-gene MLP | VAE (encoder-decoder) | cINN + VAE | Deep generative + Perturb-adapter | VAE-GAN + GCN | Neural OT (ICNN) |
| **핵심 혁신** | Gene regulatory network를 GNN에 통합 | Perturbation 효과의 가산적 분해 (disentangle) | Perturbation→cell state 양방향 매핑 | FCFP fingerprint×dose 임베딩 | GCN으로 sparse/noisy 데이터 보강 | 최적 수송 이론으로 분포 매핑 |
| **Perturbation 인코딩** | Gene embedding + perturbation embedding (GNN node) | Latent perturbation embedding (dose-dependent scaling) | ChemicalVAE/GenotypeVAE → cINN conditioning | rFCFP: FCFP×log10(dose+1) | VAE latent (disentangled) | Perturbation별 개별 transport map |
| **약물 처리** | ❌ 없음 | ✅ chemCPA 확장 (chemical descriptor) | ✅ ChemicalVAE (SMILES→latent, ZINC 250K 사전학습) | ✅ rFCFP (SMILES→fingerprint) | ✅ CMAP 기반 간접 처리 | ✅ 34개 화학 물질 평가 |
| **조합 perturbation** | ✅ 2-gene 조합 예측 | ✅ 가산적 분해로 조합 처리 | ✅ GenotypeVAE (GO annotation, ~177M 조합 학습) | ❌ 단일 perturbation | ❌ 단일 perturbation | ❌ 단일 perturbation |
| **지식 그래프** | ✅ Gene coexpression + GO | ❌ 없음 | ✅ GO annotations (15,988 terms, 18,832 genes) | ❌ 없음 | ✅ Gene regulatory network (iDREM) | ❌ 없음 |
| **Dose-response** | ❌ 없음 | ✅ Dose-dependent scaling | ✅ ChemicalVAE가 dose 처리 | ✅ log10(dose+1) 곱셈 | ❌ 없음 | ❌ 없음 |

### 1.2 데이터셋 및 정량 결과 비교

| 차원 | GEARS | CPA | PerturbNet | PRnet | UNAGI | CellOT |
|------|-------|-----|------------|-------|-------|--------|
| **주요 학습 데이터** | Replogle 2022, Norman 2019 | sci-Plex, Norman, Thompson | sci-Plex, LINCS-Drug, Norman, coding variants | L1000 bulk (~100M), sci-Plex (188 compounds) | IPF lung (54 samples, ~231K cells), COVID | Melanoma 4i (34 drugs), scRNA-seq (9 treatments) |
| **학습 데이터 규모** | ~2.5M cells (Replogle) | ~650K cells (sci-Plex) | sci-Plex 648K cells, LINCS 1.3M cells | ~100M bulk + 648K single-cell | ~231K cells | ~수만 cells |
| **sci-Plex R²** | — | 0.81 (drug combinations) | **0.984** (median) | **0.969** (unseen compounds) | — | — |
| **LINCS-Drug R²** | — | 0.857–0.891 (chemCPA) | **0.874–0.933** (median) | PCC=0.8 (unseen compounds) | — | — |
| **Norman R²** | 0.858 (all genes) | — | **0.942** (all genes) | — | — | — |
| **Unseen perturbation** | 2x Pearson over baselines | 97.6% missing combinations imputed | ✅ LINCS/학습되지 않은 화합물 예측 | ✅ PCC=0.8 unseen compounds | ✅ CMAP 기반 zero-shot | ✅ Unseen patients, cross-species |
| **MoA 직접 예측** | ❌ | ❌ (단, drug embedding이 MoA별 클러스터링) | ❌ | ✅ GSEA-based drug screening for 233 diseases | ✅ CMAP 기반 disease→drug | ❌ (단, transport cost로 MoA 계층적 군집화) |
| **실험 검증** | ❌ | ❌ | ❌ | ✅ SCLC IC50 < 10 μM, CRC natural compounds | ✅ Nifedipine anti-fibrotic PCLS | ✅ MEKi 저항성 메커니즘 발견 |

### 1.3 한계점 비교

| 모델 | 핵심 한계 | MoA 예측 관련 한계 |
|------|----------|-------------------|
| **GEARS** | Knowledge graph에 poorly connected gene 예측 불가; drug perturbation 미지원 | 약물 perturbation 자체를 처리할 수 없어 MoA 예측 불가 |
| **CPA** | Latent space 해석이 제한적; 약물 embedding 클러스터링은 MoA 분류와 동일하지 않음 | MoA 클러스터링 발견은 있으나, 분류/예측 task로 발전시킨 연구 없음 |
| **PerturbNet** | Unseen perturbation on unseen cell type 예측 불가; 학습에 대량 GO annotation 필요 | Drug MoA 분류 task 수행 안 함; chemical+genetic 동시 예측은 하지만 cross-perturbation matching 없음 |
| **PRnet** | Bulk L1000 기반 주 학습; single-cell은 보조; GSEA 해석의 생물학적 타당성 검증 제한적 | 가장 MoA에 근접하나, GSEA 기반이라 single-cell 해상도의 이점 활용 부족; CRISPR perturbation과의 matching 없음 |
| **UNAGI** | CMAP (bulk level) 기반 drug screening; unsupervised라 perturbation ground truth 불필요하지만 검증 어려움 | Perturb-seq CRISPR 데이터 미활용; CMAP bulk 한계; disease-driven이지 MoA 분류가 아님 |
| **CellOT** | Per-perturbation 개별 모델 (공유 불가); 강한 perturbation 시 실패; 결정론적 trajectory만 | MoA 예측 task 없음; perturbation 간 공유 메커니즘 학습 불가 (개별 모델); 확장성 제한 |

---

## 2. 핵심 데이터셋 분석

### 2.1 데이터셋 상세 명세

| 데이터셋 | Perturbation 유형 | 규모 | Cell line | MoA 관련성 | 접근성 |
|----------|-------------------|------|-----------|------------|--------|
| **Replogle 2022** | CRISPRi (genome-scale) | >2.5M cells, 5,000+ genes | K562 (주) | ⭐⭐⭐ 전 유전자 KO profile = MoA 매칭 기반 | GEO (추정) |
| **Norman 2019** | CRISPRa/CRISPRi (2-gene 조합) | ~218K cells, 131 조합 | K562 | ⭐⭐ 조합 효과 → pathway interaction | GEO 확인 필요 |
| **sci-Plex** | 약물 (188 compounds, 4 doses) | ~650K cells | 3 cancer cell lines (A549, K562, MCF7) | ⭐⭐⭐⭐⭐ 직접적 drug perturbation data | GEO (PerturbNet에서 648,737 cells로 전처리) |
| **LINCS L1000** | 약물 (수만 compounds) | ~100M bulk observations | 88 cell lines | ⭐⭐⭐⭐ 대규모 drug signature (bulk) | 공개 (PRnet에서 175,549 compounds 활용) |

### 2.2 sci-Plex 상세 (타 논문에서 인용된 수치 기반)

> ⚠ 전문 파싱 불가 — PerturbNet/PRnet 논문에서 인용된 전처리 수치 기반

- **원시 규모**: 188 compounds × 3 cell lines × 4 doses ≈ ~2,256 조건
- **PerturbNet 전처리 후**: 648,737 cells, 5,087 genes, 180 drug treatments (filtering 후)
- **PRnet 전처리 후**: 188 compounds 유지, single-cell 해상도
- **Cell lines**: A549 (폐암), K562 (백혈병), MCF7 (유방암)
- **주요 약물군**: HDAC inhibitor, CDK inhibitor, proteasome inhibitor 등
- **MoA annotation**: DrugBank/ATC code 기반 분류 가능

### 2.3 데이터 접근성 요약

- **Replogle 2022**: GEO accession 확인 필요 (PerturbNet 논문에서 참조). 전문 파싱 불가로 정확한 accession 미확인
- **Norman 2019**: GEO accession 확인 필요. GEARS 논문에서 참조
- **sci-Plex**: PMID 31806696, PerturbNet/CPA/PRnet에서 공통 벤치마크로 사용 → 전처리 코드 공개 가능성 높음
- **LINCS L1000**: 공개 접근 가능 (clue.io), PRnet이 175,549 compounds로 활용

---

## 3. 연구 갭 정밀 분석

### 3.1 확인된 연구 갭

> **핵심 갭**: CRISPR perturbation profile과 drug perturbation profile을 single-cell 해상도에서 직접 매칭하여 MoA를 예측하는 연구가 **존재하지 않음**

| 연구 방향 | 기존 연구 상태 | 갭 |
|-----------|---------------|-----|
| CRISPR perturbation 효과 예측 | GEARS, PerturbNet (genetic) | — |
| Drug perturbation 효과 예측 | CPA, PerturbNet (chemical), PRnet | — |
| CRISPR ↔ Drug cross-perturbation | ❌ 없음 | **이 연구의 타겟** |
| Drug MoA 분류/예측 | PRnet (GSEA), UNAGI (CMAP) | Bulk 수준에 머물름 |
| Single-cell 해상도 MoA | ❌ 없음 | **이 연구의 타겟** |
| Perturbation latent space에서 MoA 클러스터링 | CPA (발견적), CellOT (계층적) | 분류 task로 발전시킨 연구 없음 |

### 3.2 갭의 타당성 검증

1. **CPA의 drug embedding 클러스터링** [전문 파싱 근거]: CPA 논문에서 약물 embedding이 MoA별로 클러스터링됨을 보였으나, 이를 MoA 분류 task로 활용한 후속 연구 없음 → **latent space의 MoA 정보 활용 미개척**

2. **PerturbNet의 chemical+genetic 동시 처리** [전문 파싱 근거]: ChemicalVAE와 GenotypeVAE가 별도 학습되어, 동일 latent space에서 CRISPR-약물 matching이 직접 수행되지 않음 → **cross-perturbation latent space 미구축**

3. **PRnet의 GSEA 한계** [전문 파싱 근거]: GSEA는 pathway enrichment 기반으로, single-cell의 세포 간 이질성을 반영하지 못함 → **single-cell 해상도의 이점 미활용**

4. **CellOT의 transport cost clustering** [전문 파싱 근거]: 약물 간 transport cost의 계층적 군집화가 MoA와 유사한 그룹을 형성하지만, (a) per-perturbation 모델로 확장성 제한, (b) MoA 분류 task로 정식화되지 않음

### 3.3 갭의 중요성

- **실제 수요**: 신약 개발에서 MoA 규명은 필수적. 임상 시험 실패의 주요 원인 중 하나가 MoA 불명확
- **기술적 가능성**: PerturbNet R²=0.984 (sci-Plex), PRnet PCC=0.8 (unseen compounds) 등, perturbation 효과 예측이 이미 고정확도 → latent space 품질 충분
- **단일 세포의 추가 가치**: Bulk CMap은 약물 반응의 이질성(예: 저항성 subpopulation)을 포착하지 못함. Single-cell 해상도에서는 subpopulation-specific MoA 파악 가능

---

## 4. 제안 방법론 검토: 기존 모델 시사점 기반 구체화

### 4.1 PerturbNet의 시사점 → 아키텍처 기반

**채택 가능 요소**:
- **cINN (conditional Invertible Neural Network)**: Perturbation space ↔ cell state space의 양방향 매핑. 기존 VAE의 한 방향 병목 해결
- **GenotypeVAE**: GO annotation 기반 유전자 기능 임베딩 (~177M 조합 학습). CRISPR perturbation을 기능적 공간에 배치
- **ChemicalVAE**: SMILES → latent. 약물을 화학 구조 공간에 배치

**제안 수정**: GenotypeVAE + ChemicalVAE를 **공통 latent space**로 통합 (PerturbNet은 별도 학습). 이 공간에서 CRISPR KO와 drug perturbation이 같은 좌표계를 가져야 MoA 매칭 가능.

### 4.2 CPA의 시사점 → Disentanglement 기반

**채택 가능 요소**:
- **가산적 perturbation 분해**: basal state + perturbation effect + covariate → perturbation 효과만 분리 가능
- **Adversarial training**: perturbation embedding이 covariate 정보를 포함하지 않도록 보장
- **Drug embedding의 MoA 클러스터링**: 이미 검증된 현상

**제안 수정**: CPA의 disentangle 구조를 채택하되, perturbation embedding에 **MoA-aware contrastive loss**를 추가. 같은 MoA의 약물/CRISPR perturbation은 가깝게, 다른 MoA는 멀게.

### 4.3 PRnet의 시사점 → Drug 임베딩 + 실증 검증 기반

**채택 가능 요소**:
- **rFCFP 임베딩**: FCFP fingerprint × log10(dose+1) → 화학 구조 + 용량 정보를 간결하게 인코딩
- **GSEA-based drug screening**: 233개 질환에 대한 drug atlas 구축 → MoA 예측 검증 프레임워크로 활용 가능
- **실험 검증 파이프라인**: IC50 측정으로 in silico 예측의 생물학적 타당성 검증

**제안 수정**: rFCFP를 drug encoder로 채택. 단, GSEA 대신 **latent space similarity**로 대체 (single-cell 해상도 유지). PRnet의 실험 검증 결과는 baseline 비교에 활용.

### 4.4 CellOT의 시사점 → 이질성 포착 기반

**채택 가능 요소**:
- **최적 수송 기반 분포 매핑**: 단순 mean shift가 아닌 full distribution mapping
- **Subpopulation-specific 약물 반응**: 저항성 subpopulation 포착

**한계 및 대안**: Per-perturbation 모델은 확장성 제한. 대신 **conditional OT** (perturbation embedding으로 condition)를 사용하여 shared model 구축.

### 4.5 통합 제안 아키텍처

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
│  │ GenotypeVAE  │         │ ChemicalVAE  │           │
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
│  Loss: Reconstruction + Adversarial + Contrastive    │
│        (MoA-aware: same MoA → close, diff → far)     │
└─────────────────────────────────────────────────────┘
```

**핵심 설계 결정**:

1. **GenotypeVAE + ChemicalVAE → Shared Latent Space**: PerturbNet의 별도 학습과 달리, 공통 latent space에서 CRISPR-약물 직접 매칭
2. **rFCFP 임베딩 (PRnet)**: Drug encoder로 chemical fingerprint × dose 채택. ChemicalVAE (SMILES→latent)보다 해석 가능성 높음
3. **CPA-style disentanglement**: Basal state와 perturbation effect 분리. Adversarial loss로 perturbation embedding이 covariate-independent 보장
4. **MoA-aware contrastive loss**: [새로 제안] 같은 MoA 그룹의 약물/CRISPR perturbation을 latent space에서 가깝게. Hard negative mining: 같은 pathway 내 다른 MoA
5. **3개 출력 헤드**:
   - MoA Classification Head: 약물의 MoA 라벨 분류 (ATC code 기반)
   - Cell State Decoder: Perturbation 효과 재구성 (평가용)
   - Cross-Perturbation Matching: CRISPR KO ↔ Drug profile 유사도 (MoA 추론)

---

## 5. 평가 프로토콜 구체화

### 5.1 데이터 분할 전략

| 전략 | 설명 | 평가 목적 |
|------|------|----------|
| **Leave-compound-out** | 특정 약물을 test에 배제 | Unseen compound에 대한 MoA 예측 |
| **Leave-MoA-out** | 특정 MoA 클래스 전체를 test에 배제 | Zero-shot MoA 분류 |
| **Leave-cell-line-out** | 특정 cell line을 test에 배제 | Cross-cell-line 일반화 |
| **Cross-perturbation** | CRISPR KO로 학습 → 약물 MoA 예측 | 핵심: CRISPR→Drug 전이 |

### 5.2 평가 지표

| 지표 | 용도 | Baseline |
|------|------|----------|
| **MoA 분류 정확도** (Top-1/3/5) | MoA 예측 성능 | Random, CMap (bulk), PRnet (GSEA) |
| **R² (gene-level)** | Perturbation 효과 재구성 | GEARS, CPA, PerturbNet 원본 수치 |
| **Pearson correlation** | Latent space 유사도 타당성 | PRnet PCC=0.8 기준 |
| **Adjusted Rand Index (ARI)** | Latent space MoA 클러스터링 품질 | CPA 클러스터링 기준 |
| **AUROC (MoA binary)** | 특정 MoA 탐지 | CMap connectivity score |

### 5.3 Baseline 비교 대상

| Baseline | 비교 포인트 | 예상 결과 |
|----------|------------|----------|
| **CMap (bulk)** | Bulk 대비 single-cell 이점 | 본 제안이 subpopulation-specific MoA에서 우세 |
| **CPA (원본)** | Drug embedding 클러스터링 | 본 제안이 contrastive loss로 MoA 분류에서 우세 |
| **PerturbNet** | Perturbation 예측 정확도 | R²은 PerturbNet이 높을 수 있으나, MoA 분류는 본 제안만 수행 |
| **PRnet (GSEA)** | Drug screening 접근 비교 | GSEA는 bulk 기반; 본 제안은 single-cell heterogeneity 반영 |
| **Random Forest (gene expression → MoA)** | 단순 분류기 baseline | Feature engineering 없이 raw expression → MoA의 upper bound 확인 |

---

## 6. 리스크 업데이트 (전문 분석 기반)

| 리스크 | 확률 | 심각도 | 완화 전략 | 근거 |
|--------|------|--------|----------|------|
| **CRISPR KO ≠ Drug inhibition** | **높음** | 높음 | CRISPRi (partial knockdown) 데이터 활용; stoichiometry 보정 항 추가 | GEARS/CPA 모두 이 문제를 직접 다루지 않음. Replogle 2022가 CRISPRi 사용 → KO보다 drug inhibition에 더 근접 |
| **Latent space 정렬 불가** | 중간 | 높음 | Warm-start: GenotypeVAE/ChemicalVAE를 각각 사전학습 후 공통 space로 미세조정 | PerturbNet은 별도 VAE 학습 → 공통 space가 가능한지 미검증 |
| **sci-Plex 3 cell line 일반화 부족** | 중간 | 중간 | LINCS L1000 (88 cell lines, bulk) 보조 학습; cross-cell-line 평가 추가 | PRnet이 88 cell lines 활용했으나 bulk 수준 |
| **MoA 라벨 불균형** | 중간 | 중간 | Focal loss; few-shot 대비 metric learning | sci-Plex 188 compounds의 MoA 분포 확인 필요 |
| **Contrastive loss 해킹** (trivial solution) | 낮음–중간 | 중간 | Temperature scaling; hard negative mining; alignment+uniformity loss | SimCLR/MoCo 문헌에서 검증된 완화법 |
| **Novelty 확보** | 중간 | 높음 | Cross-perturbation contrastive learning이 핵심 차별화; CRISPR→Drug zero-shot MoA는 선행 연구 없음 | 본 분석에서 확인된 갭 |

---

## 7. 즉각적 다음 단계 (우선순위순)

1. **sci-Plex 데이터 전처리 파이프라인 구축**
   - PerturbNet/CPA/PRnet의 공개 코드에서 sci-Plex 전처리 재현
   - MoA annotation (DrugBank/ATC code) 매핑
   - 근거: 3개 모델이 공통으로 사용 → 전처리 코드 공개 가능성 높음

2. **CPA drug embedding MoA 클러스터링 재현**
   - CPA 공식 코드로 sci-Plex 학습 → drug embedding t-SNE/UMAP
   - MoA별 클러스터링 정도 정량화 (ARI, NMI)
   - 근거: CPA 논문이 클러스터링 보고했으나 정량화 없음 → 이 연구의 baseline

3. **Cross-perturbation contrastive learning 프로토타입**
   - PyTorch + PyG 구현
   - GenotypeVAE + rFCFP drug encoder → shared latent space
   - MoA-aware contrastive loss 실험

4. **Replogle 2022 데이터 접근**
   - GEO accession 확인 후 다운로드 테스트
   - CRISPRi (KO가 아님) 데이터 = drug inhibition과 더 유사 → 리스크 완화에 기여

---

## 8. 모델별 상세 분석 요약

### GEARS
- **근거 수준**: 전문 파싱 완료
- **핵심 기여**: Gene regulatory network (coexpression + GO)를 GNN에 통합하여 perturbation 효과 예측. Cross-gene MLP가 gene 간 dependency 포착.
- **MoA 관련성**: 낮음. Drug perturbation 미지원. 단, GNN 구조와 gene embedding 방식은 참고 가능.
- **한계**: Knowledge graph에 poorly connected gene 예측 불가; drug 미지원.

### CPA
- **근거 수준**: 전문 파싱 완료
- **핵심 기여**: Perturbation 효과의 가산적 분해 (basal + perturbation + covariate). Adversarial training으로 disentanglement 보장. chemCPA로 chemical descriptor 도입.
- **MoA 관련성**: 중간. Drug embedding이 MoA별 클러스터링됨을 보였으나, 이를 MoA 분류 task로 활용하지 않음. **본 연구의 직접적 출발점**.
- **한계**: Latent space 해석 제한적; 클러스터링 → 분류 발전 필요.

### PerturbNet
- **근거 수준**: 전문 파싱 완료
- **핵심 기여**: cINN으로 perturbation↔cell state 양방향 매핑. 가장 높은 R² (sci-Plex 0.984). Chemical+genetic perturbation 모두 처리.
- **MoA 관련성**: 중간. Chemical+genetic을 모두 처리하지만, cross-perturbation matching 없음. 아키텍처 구성요소 (cINN, VAE)는 직접 활용 가능.
- **한계**: Unseen perturbation on unseen cell type 불가; 대량 GO annotation 필요.

### PRnet
- **근거 수준**: 전문 파싱 완료
- **핵심 기여**: rFCFP 임베딩 (화학 구조 + 용량); GSEA 기반 233 질환 drug screening; SCLC/CRC 실험 검증. Bulk L1000 ~100M 관측치 활용.
- **MoA 관련성**: 높음. 가장 직접적으로 약물 발견 목표. 단, GSEA 기반이라 single-cell 이점 미활용.
- **한계**: Bulk가 주 데이터; single-cell은 보조; CRISPR perturbation과의 연결 없음.

### UNAGI
- **근거 수준**: 전문 파싱 완료
- **핵심 기여**: VAE-GAN+GCN으로 sparse/noisy scRNA-seq 보강; 질병 진행 궤적 + 약물 스크리닝 통합; unsupervised 접근.
- **MoA 관련성**: 낮음–중간. Disease-driven drug screening이지 MoA 분류가 아님. CMAP (bulk) 기반.
- **한계**: Perturb-seq CRISPR 데이터 미활용; CMAP bulk 한계; unsupervised 검증 어려움.

### CellOT
- **근거 수준**: 전문 파싱 완료
- **핵심 기여**: Neural OT로 perturbation 전후 분포 매핑. Subpopulation-specific 반응 포착. Cross-species/patient 일반화.
- **MoA 관련성**: 낮음. Transport cost로 MoA 유사 군집화 가능하나, MoA 분류 task 없음. Per-perturbation 모델로 확장성 제한.
- **한계**: Per-perturbation 개별 모델; 강한 perturbation 시 실패; 결정론적만.

---

*이 보고서는 6편의 전문 PDF 파싱 결과와 3편의 초록/메타데이터를 기반으로 작성되었습니다. 구체적 수치는 해당 논문의 실험 결과에서 직접 인용하였으며, [미검증] 또는 [초록 수준] 표시가 없는 항목은 전문 파싱으로 확인된 근거입니다.*
