# Stage 03 — Planning

## Workflow
1. `docs/03_planning.md` 가이드 확인
2. `stages/02_framing.md`에서 베이스라인/타겟 확인
3. 실험 설계 작성 (모델, 피처, 하이퍼파라미터, 평가 전략)
4. **사용자 컨펌 획득** (타겟 성능 + 실험 설계)
5. 산출물 → `outputs/planning/run_XX/`에 저장
6. 아래 지식 업데이트 (과거 run에서 검증된 인사이트만 통합)
7. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 실험 설계서 작성 완료
- 사용자 컨펌 획득 (`confirmed by user [날짜]`)

---

## 연구 질문 (stages/02에서)

**"Factorized causal representations를 Independent Causal Mechanism 원리로 정규화했을 때, 섭동 효과가 세포 유형에 걸쳐 불변하는 인과 모듈로 분해되며, 이 모듈들의 조합이 타겟 세포 유형의 학습 데이터 없이 조합 섭동 효과를 예측할 수 있는가?"**

---

## 1. 방법론 아키텍처

```
Module 1: FCR 인코더 (z_x, z_t, z_tx 분해)
  → Module 2: ICM 정규화 (z_tx 불변성 강제)
    → Module 3: 모듈 분해 (z_tx → pathway-level modules)
      → Module 4: 조합 함수 (모듈 조합 → 조합 섭동 예측)
        → Module 5: 교차 세포 유형 전이 (불변 모듈 → 제로샷)
```

의존: 1→2→3→4→5 순차. 각 모듈은 소거 실험에서 개별 기여도 측정 가능.

## 2. 설계 결정과 근거

| 설계 결정 | 근거 | 대안 검토 |
|-----------|------|----------|
| FCR을 기반 표현 학습으로 채택 | z_x/z_t/z_tx 분해 + 식별 가능성 증명(arXiv:2410.22472) | CPA(조합적이나 분해 없음), SAMS-VAE(가법 분해 한계) |
| ICM 정규화로 불변성 강제 | Schölkopf group IEM(arXiv:2406.14302) — 메커니즘 자율성 원리 | IRM(최적화 어려움), DRO(보수적) |
| 경로 수준 모듈 분해 | scBIG(arXiv:2602.04901)이 유전자 프로그램 구조의 효용 입증 | Gene-level(너무 미세), Cell-level(너무 거시) |
| Norman + Replogle 데이터 | 조합 섭동(Norman) + 다세포 유형(Replogle) 동시 확보 | sci-Plex(약물 섭동, CRISPR와 다른 메커니즘) |
| R2 + DEG + ARI 복합 평가 | 단일 지표로는 불충분(Problem 1: 평가 지표 붕괴) | MSE만 사용(Shesha가 보여준 한계) |

## 3. 실험 설계

### Module 1: FCR 인코더
- **아키텍처**: 인코더 q(z|x,t) → z_x(공변량), z_t(처리), z_tx(상호작용) 분해
- **학습**: FCR 원논문 구현 기반, block-wise identifiability 보장 조건 충족
- **입력**: 단일 세포 발현 벡터 x ∈ R^G + 섭동 식별자 t
- **하이퍼파라미터**: z_dim=[16,32,64], β-VAE weight=[0.1,1.0,10.0]

### Module 2: ICM 정규화
- **목적**: z_tx가 세포 유형에 걸쳐 불변하도록 정규화
- **정규화 항**: L_ICM = Σ_c ||∇_c f_tx(z_tx; c)||² (맥락 c에 대한 기울기 패널티)
  - c = 세포 유형 원-핫 인코딩
  - f_tx = z_tx에서 섭동 효과를 디코딩하는 함수
- **대안**: MMD(z_tx^A, z_tx^B) — 분포 정렬, 더 약한 제약
- **소거 실험**: ICM 정규화 유무별 z_tx 상관관계 비교

### Module 3: 모듈 분해
- **z_tx → {m_1, ..., m_K}**: K개 경로 수준 모듈로 희소 분해
- **구현**: Sparse attention 또는 hard assignment (Gumbel-Softmax)
- **경로 사전**: MSigDB Hallmark(50경로) + KEGG(~300경로)로 초기화, 미세조정 허용
- **평가**: 학습된 모듈 ↔ 알려진 경로의 ARI

### Module 4: 조합 함수
- **h(m_A, m_B) → Δy_{AB}**: 두 모듈의 조합 효과 예측
- **구현 옵션**:
  - (a) 가법: m_A + m_B (하한 baseline)
  - (b) 곱법: m_A ⊙ m_B (경로 내 상호작용)
  - (c) 학습 가능: MLP(m_A, m_B) (데이터 기반)
- **학습**: 단일 섭동 m_i로부터 조합 규칙 학습
  - 기대: 독립 경로 = 가법, 동일 경로 = 곱법/학습
- **평가**: Norman double-KO 데이터로 조합 예측 R2

### Module 5: 교차 세포 유형 전이
- **소스**: K562 (섭동 데이터 풍부)
- **타겟**: RPE1 (섭동 데이터 없이 예측)
- **방법**: 소스에서 학습된 불변 모듈 + 조합 함수를 타겟에 직접 적용
- **z_x만 타겟 세포 유형에서 학습** (공변량 분포는 다르므로)
- **평가**: RPE1 제로샷 예측 R2, z_tx 불변성 검증(r > 0.7)

## 4. 평가 설계

| 실험 | RQ | 비교 대상 | 지표 | 타겟 |
|------|-----|----------|------|------|
| FCR 분해 품질 | RQ1 | FCR(no ICM) vs FCR+ICM | z_tx cross-cell correlation | r > 0.7 (ICM) vs r < 0.3 (no ICM) |
| 조합 예측 (독립 경로) | RQ2 | Additive vs h(m_A,m_B) | R2 on double-KO | R2 > 0.5 |
| 조합 예측 (동일 경로) | RQ2 | Additive vs h(m_A,m_B) | R2 on double-KO | R2 > 0.3 (상위성 존재) |
| 제로샷 전이 | RQ3 | Scratch vs transfer | R2 on RPE1 | R2 > 0.3 |
| 모듈-경로 정렬 | — | Random vs learned | ARI vs MSigDB | ARI > 0.4 |
| 상위성 탐지 | RQ2 | 조합 잔차 | Epistasis recall | > 0.6 |

## 5. 소거 실험 매트릭스

| 구성 | Module 1 | Module 2 | Module 3 | Module 4 | 예상 결과 |
|------|----------|----------|----------|----------|----------|
| Full model | FCR | ICM | Pathway modules | Learned h | 최고 성능 |
| w/o ICM | FCR | — | Pathway modules | Learned h | 교차 세포 유형 불변성 상실 |
| w/o modules | FCR | ICM | — (z_tx 직접 사용) | Learned h | 조합 예측 저하 |
| w/o composition | FCR | ICM | Pathway modules | Additive only | 동일 경로 조합 실패 |
| CPA baseline | CPA | — | — | — | OOD 급락 |
| GEARS baseline | GEARS | — | — | — | 교차 세포 유형 미지원 |

## 6. 논문 구조 제안

1. **Introduction**: Perturb-seq의 두 최대 격차 — 조합 폭발 + 교차 세포 유형 전이 — 를 단일 프레임워크로 해결
2. **Related Work**: 섭동 예측(GEARS, CPA, CellOT), 인과 표현 학습(FCR, IEM), 조합 일반화(scBIG, MapPFN)
3. **Methods**: FCR+ICM 인코더 → 경로 모듈 분해 → 학습 가능 조합 함수 → 제로샷 전이
4. **Results**:
   - 4.1 RQ1: ICM이 z_tx 불변성을 보장하는가?
   - 4.2 RQ2: 모듈 조합이 조합 섭동을 예측하는가?
   - 4.3 RQ3: 불변 모듈이 제로샷 전이를 가능하게 하는가?
   - 4.4 상위성 탐지: 조합 잔차의 생물학적 해석
5. **Discussion**: 인과 불변성과 조합성의 시너지, 한계(식별 가능성 가정), 확장(약물 섭동)

## 7. 핵심 리스크와 완화

1. **FCR 식별 가능성 가정이 현실 데이터에서 성립 불가** → Block-wise identifiability는 약한 조건; 성립하지 않으면 z_tx가 진짜 인과 표현이 아닐 수 있음. 완화: 소거 실험으로 ICM 정규화 효과를 경험적으로 검증
2. **Norman 데이터셋 규모 제약** (~100K cells, 104 double KO) → 교차 검증으로 최대 활용. Replogle은 더 크나 조합 섭동 없음
3. **IRM-style 패널티의 최적화 어려움** → MMD 정규화를 약한 대안으로 준비. 둘 다 실패하면 단순 fine-tuning baseline
4. **상위성이 비조합적일 수 있음** → 부분적 조합 구조조차 가치. 실패 모드 자체가 상위성 탐지 도구로 기능

## 8. 데이터 확보 계획

| 데이터셋 | 접근성 | 전처리 | 비고 |
|----------|--------|--------|------|
| Norman et al. (2019) | GEO GSE133344 | Scanpy 표준 파이프라인 | 조합성 검증 핵심 |
| Replogle et al. (2022) | GEO GSE142398 | K562 + RPE1 분리 | 교차 세포 유형 전이 |
| MSigDB Hallmark | 공개 | 경로 gene set | 모듈 평가 기준 |
| KEGG/Reactome | 공개 | 경로 토폴로지 | 모듈 초기화 |

## 9. 컴퓨팅 자원 추정

| 단계 | 예상 시간 | 비고 |
|------|----------|------|
| 데이터 전처리 (Norman+Replogle) | 2-4h | Scanpy, 메모리 16GB+ |
| FCR 인코더 학습 | 1-3h/실험 | GPU 권장, RTX 4060 Ti |
| ICM 정규화 실험 | 2-4h/실험 | 추가 loss term만 |
| 조합 함수 학습 | 1-2h | 단일 KO → double KO |
| 제로샷 전이 평가 | 1h | inference only |
| 전체 소거 실험 | 1-2일 | 6개 구성 × 3-fold CV |

---

## Run 이력 (세부 내용은 outputs/planning/run_XX/ 참조)
- **run_01** (2026-04-26): NAP 기반 설계. Phase 3-4 기술실현성 심층 분석
- **run_02** (2026-04-27): NAP novelty 부족으로 루프백. 인과 조합 가설 기반 전면 재설계
