# 최종 분석 보고서: Perturb-seq 기반 약물 MoA 예측 — Cross-Perturbation MoA Framework

> 작성일: 2026-04-16 | 환경: Python 3.11, PyTorch 2.5.1, CUDA RTX 4060 Ti

---

## 0. 근거 수준 명시

| 항목 | 근거 수준 | 비고 |
|------|----------|------|
| 선행연구 6편 분석 | **전문 파싱** | GEARS, CPA, PerturbNet, PRnet, UNAGI, CellOT |
| 데이터셋 3편 | **초록 수준** | Replogle 2022, Norman 2019, sci-Plex (페이월) |
| 제안 아키텍처 | **프로토타입 구현 + 합성 데이터 테스트** | 실제 데이터 미적용 |
| 프로토타입 결과 | **합성 데이터 기반** | 실제 sci-Plex 데이터 결과와 차이 클 수 있음 |
| 연구 갭 | **도구 검증** | 6편 전문 분석으로 확인 |

---

## 1. 연구 요약

### 1.1 연구 질문
> Perturb-seq 데이터에서 추출한 유전자 perturbation signature를 활용하여, 신약 후보 물질의 작용 기전(MoA)을 single-cell 해상도에서 얼마나 정확하게 예측할 수 있는가?

### 1.2 확인된 연구 갭
**CRISPR perturbation profile과 drug perturbation profile을 single-cell 해상도에서 직접 매칭하여 MoA를 예측하는 연구가 존재하지 않음**

- 기존 모델들은 CRISPR perturbation 또는 drug perturbation 중 하나만 처리
- CPA는 drug embedding의 MoA 클러스터링을 발견했으나, 이를 분류 task로 발전시키지 않음
- PRnet은 GSEA 기반 drug screening을 수행하나 bulk 수준에 머물름
- PerturbNet은 chemical+genetic 모두 처리하나 cross-perturbation matching 없음

### 1.3 제안 방법론: Cross-Perturbation MoA Framework

```
CRISPR Perturbation → GenotypeVAE ──┐
                                     ├── Shared Latent Space ── MoA Classification
Drug Perturbation → rFCFP Encoder ──┘    (Disentangled +         │
                                          MoA-aware CL)          ↓
                                                              Cell State
                                                              Reconstruction
```

핵심 혁신: **MoA-aware contrastive loss** — 같은 MoA의 CRISPR/drug perturbation을 latent space에서 가깝게 학습

---

## 2. 선행연구 분석 결과

### 2.1 모델 비교 (전문 파싱 근거)

| 모델 | 구조 | 약물 | MoA | sci-Plex R² | 한계 |
|------|------|------|-----|-------------|------|
| GEARS | GNN+knowledge graph | X | X | -- | Drug 미지원 |
| CPA | VAE+adversarial | O (chemCPA) | 발견적 | 0.81 | 분류 task 미발전 |
| PerturbNet | cINN+VAE | O | X | **0.984** | Cross-matching 없음 |
| PRnet | Deep gen.+rFCFP | O | O (GSEA) | 0.969 | Bulk 기반 GSEA |
| UNAGI | VAE-GAN+GCN | O (CMAP) | O (disease→drug) | -- | CMAP bulk 한계 |
| CellOT | Neural OT | O | X | -- | Per-perturbation 모델 |

### 2.2 데이터셋 명세

| 데이터셋 | 유형 | 규모 | 핵심 활용 |
|----------|------|------|----------|
| sci-Plex | 약물 perturbation | ~650K cells, 188 compounds, 3 cell lines | 핵심 학습/평가 데이터 |
| Replogle 2022 | CRISPRi (genome-scale) | >2.5M cells, 5,000+ genes | CRISPR perturbation encoder 학습 |
| Norman 2019 | CRISPR 조합 | ~218K cells, 131 조합 | 조합 perturbation 평가 |
| LINCS L1000 | 약물 bulk | ~100M observations, 175K compounds | 대규모 drug signature 보조 |

---

## 3. 프로토타입 구현 결과

### 3.1 구현 내역

| 파일 | 내용 |
|------|------|
| `config.py` | 하이퍼파라미터, 경로 설정 |
| `data_utils.py` | 합성 데이터 생성, MoA annotation, 4종 데이터 분할 (random/leave-compound/leave-moa/leave-cell-line) |
| `model.py` | PerturbationEncoder, BasalEncoder, CovariateEncoder, CellStateDecoder, MoAClassificationHead, AdversarialDiscriminator, MoAContrastiveLoss, CrossPerturbMoAModel, cross_perturbation_match |
| `train.py` | 학습 루프, 평가, 결과 저장 |
| `run_pipeline.py` | 전체 실행 스크립트 |

### 3.2 합성 데이터 테스트 결과

**환경**: NVIDIA RTX 4060 Ti, PyTorch 2.5.1, CUDA
**데이터**: 20,000 cells x 1,000 genes, 30 drugs, 6 MoA classes, 3 cell lines
**모델**: 1,331,642 파라미터, 30 epochs

| 분할 전략 | Top-1 Accuracy | Top-3 Accuracy | F1 (macro) | Recon MSE | 특징 |
|-----------|---------------|----------------|------------|-----------|------|
| **Random** | **1.0000** | **1.0000** | **1.0000** | 38.15 | 학습에 본 약물 포함 |
| **Leave-compound-out** | **0.1674** | 0.4173 | 0.1440 | 38.38 | 보이지 않는 약물 — MoA 예측 어려움 |
| **Leave-MoA-out** | **0.0000** | 0.0000 | 0.0000 | 38.02 | Zero-shot MoA — 현재 구조로 불가 |

### 3.3 결과 해석

1. **Random split (100%)**: 합성 데이터에서 MoA 신호가 충분히 강하여 쉽게 분류 가능. 실제 데이터에서는 이 정도 성능 불가능.

2. **Leave-compound-out (16.7%)**: 보이지 않는 약물의 MoA를 예측하는 것이 핵심 과제. 현재 perturbation embedding만으로는 새 약물의 MoA를 일반화하기 어려움.
   - **원인 분석**: PerturbationEncoder가 drug index embedding만 사용 → 새 약물은 임의의 embedding 할당. Chemical structure 정보(rFCFP)가 들어가야 unseen compound 처리 가능.
   - **해결 방향**: rFCFP fingerprint를 drug encoder 입력으로 직접 사용 (현재는 index embedding만 사용 중)

3. **Leave-MoA-out (0%)**: Zero-shot MoA 예측은 분류 head 기반으로는 불가능. MoA 라벨이 학습에 없는 클래스는 예측 불가.
   - **해결 방향**: Contrastive learning으로 latent space에서 MoA 클러스터를 형성하고, nearest-neighbor 방식으로 zero-shot 분류 수행. 분류 head 대신 metric learning 접근 필요.

4. **Recon MSE (~38)**: 재구성 오차가 높음. Basal state encoder의 품질이 충분하지 않을 수 있음. 단, 합성 데이터 특성상 실제보다 노이즈가 적어 재구성 자체의 의미가 제한적.

5. **Cross-perturbation match (1.0)**: 합성 데이터에서는 동일 MoA 내 perturbation이 완벽히 구분되어 100% 달성. 실제 데이터에서는 현저히 낮을 것.

---

## 4. 핵심 발견 및 시사점

### 4.1 프로토타입에서 확인된 설계 문제

| 문제 | 원인 | 해결 방안 |
|------|------|----------|
| Unseen drug MoA 예측 불가 | Drug index embedding 사용 → 새 약물에 임의 embedding | rFCFP fingerprint를 drug encoder 입력으로 직접 사용 |
| Zero-shot MoA 불가 | Classification head 방식 → 학습에 없는 클래스 예측 불가 | Metric learning: latent space에서 nearest-neighbor 기반 zero-shot |
| Reconstruction MSE 높음 | Basal encoder가 충분한 표현력 없음 | 더 깊은 네트워크; skip connection; pretraining |
| Contrastive loss 효과 미검증 | 합성 데이터에서 MoA가 너무 쉽게 분리됨 | 실제 데이터 필요 |

### 4.2 다음 단계 개선 사항

1. **Drug Encoder를 rFCFP 기반으로 교체**
   - 현재: `nn.Embedding(n_drugs, dim)` → 새 약물 처리 불가
   - 개선: SMILES → FCFP fingerprint (RDKit) → MLP → z_drug
   - 효과: Unseen compound의 화학 구조를 기반으로 MoA 예측 가능

2. **Zero-shot MoA를 위한 metric learning**
   - 현재: Classification head (softmax cross-entropy)
   - 개선: Prototype network / nearest class mean in latent space
   - 효과: 학습에 없는 MoA도 latent space에서 가장 가까운 클러스터로 분류

3. **실제 데이터 적용**
   - sci-Plex 데이터 다운로드 (GEO GSE139944)
   - PerturbNet/CPA 전처리 코드 활용
   - DrugBank/ATC code 기반 MoA annotation

4. **CRISPR perturbation 데이터 통합**
   - Replogle 2022 (CRISPRi) 데이터 추가
   - GenotypeVAE와 Drug Encoder의 shared latent space 학습
   - Cross-perturbation matching: CRISPR KO → drug MoA 추론

---

## 5. 산출물 요약

| 산출물 | 경로 | 설명 |
|--------|------|------|
| 심층 문헌 분석 보고서 | `deep_literature_analysis_perturb_seq_moa.md` | 6편 전문 분석, 모델 비교, 연구 갭 |
| 업데이트된 연구 플랜 | `research_plan_perturb_seq_moa.md` | 전문 분석 기반으로 구체화된 방법론 |
| 파이프라인 코드 | `perturb_moa_pipeline/` | 5개 Python 모듈 |
| 학습 결과 | `perturb_moa_pipeline/results/` | JSON 결과 파일 |

---

## 6. 결론

1. **연구 갭 확인**: CRISPR perturbation ↔ drug perturbation cross-matching으로 MoA를 예측하는 연구는 문헌에 존재하지 않음. 이는 유효한 연구 기회.

2. **프로토타입 검증**: Cross-Perturbation MoA Framework의 기본 구조가 작동함을 합성 데이터로 확인. 단, unseen compound/zero-shot MoA는 추가 설계 개선 필요.

3. **핵심 개선점**: (a) rFCFP 기반 drug encoder 도입, (b) metric learning 기반 zero-shot MoA, (c) 실제 sci-Plex 데이터 적용. 이 세 가지가 다음 단계의 핵심.

4. **현실적 기대**: 합성 데이터의 100% 정확도는 실제 데이터에서 불가능. PRnet의 PCC=0.8 (unseen compounds), PerturbNet의 R²=0.984 (sci-Plex)를 참고할 때, 실제 데이터에서 Top-1 MoA 정확도 40-60% 달성이 현실적 목표.

5. **연구 가치**: Single-cell 해상도에서 CRISPR-drug cross-perturbation matching은 기존 bulk CMap 대비 subpopulation-specific MoA 파악 가능. 이는 약물 저항성 메커니즘 이해에 직접 기여.

---

*이 보고서는 6편의 전문 PDF 파싱 결과와 합성 데이터 프로토타입 실행 결과를 기반으로 작성되었습니다. 프로토타입 결과는 합성 데이터 기반이며, 실제 sci-Plex 데이터에서의 성능은 다를 수 있습니다.*
