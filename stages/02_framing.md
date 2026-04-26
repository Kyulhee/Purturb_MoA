# Stage 02 — Framing

## Workflow
1. `docs/02_framing.md` 가이드 확인
2. `stages/01_literature_review.md`에서 인사이트 확인
3. 연구 질문, 베이스라인, 타겟 성능 수치 정의
4. 산출물 → `outputs/framing/run_XX/`에 저장
5. 아래 지식 업데이트 (검증된 인사이트만 통합)
6. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 연구 질문이 단일 문장으로 기술 가능
- 베이스라인 수치가 문헌 근거와 함께 명시됨
- 타겟 성능 수치가 베이스라인 기반으로 설정됨
- 평가 지표 및 데이터셋 확정

---

## 검증된 핵심 지식

### 연구 질문
**"대조 학습 기반 약물 임베딩이 sci-Plex 단일 세포 섭동 데이터에서 MoA 클러스터링 품질을 개선하는가?"**

### 평가 전략
| 평가 방식 | 지표 | 비고 |
|-----------|------|------|
| leave-compound-out | Top-1 accuracy, F1 macro | 분류 (보이지 않는 약물) |
| leave-MoA-out | Silhouette, ARI, NMI | **클러스터링 품질** (분류 accuracy 아님) |
| Alignment/Uniformity | alignment score, uniformity score | Wang & Isola 기반 |

### 베이스라인 계층
1. Random: ~6.25% (16개 클래스)
2. Simple DNN (GPAR-style): 978→512→256
3. PANACEA top methods
4. chemCPA: r²=0.68 (DEGs, pretrained)
5. Our method (대조 학습 + MoA-aware)

### 타겟 수치
- leave-compound-out Top-1 > 60%
- leave-MoA-out ARI > 0.30
- Alignment 50% 개선 (vs non-contrastive)

### 데이터
- sci-Plex3: GEO GSE139944 (GSM4150378)
- chemCPA 전처리 h5ad 다운로드 가능
- **16개 MoA 카테고리** (pathway_level_1, "Other" 처리 결정 필요)

### 미해결
- 실제 데이터 다운로드 + MoA 분포 정량 확인 (다운로드 진행 중이었음)
- "Other" 카테고리 제외 vs 포함 결정
- Simple DNN baseline 수치 측정 (Planning/Analysis에서)

### 다음 단계에 전달
1. Loss weight 설계: L_recon vs L_contrastive 균형 (절대값이 아닌 effective gradient 기준)
2. Drug encoder fine-tuning 허용 필수
3. 데이터 다운로드 후 MoA 분포 → 클래스 불균형 대응 전략

---

## Run 이력 (세부 내용은 outputs/framing/run_XX/ 참조)
- run_01: 연구 질문/평가 전략/베이스라인 계층/타겟 수치 정의. MoA 16개 확인. 데이터 다운로드 진행 중
