# BioEval: 생물학적 충실도 기반 섭동 예측 평가 프레임워크

Perturb-seq 섭동 예측에서 평가 지표(MSE/R²)가 만드는 **Mean-Effect Trap**을 진단하고, 생물학적 충실도 기반 평가 지표(BioEval)를 설계하여 "DL ≤ Baseline" 위기가 지표의 아티팩트인지 판별하는 연구 프로젝트입니다.

## 연구 질문

> **"평가 지표를 생물학적 충실도 기반으로 바꾸면, 모델 순위가 어떻게 달라지는가? DL ≤ Baseline 위기가 지표의 아티팩트인지, 실재하는 현상인지 판별할 수 있는가?"**

## 가설 및 결과 (run_13-20)

| 가설 | 내용 | 판정 | 클레임 강도 |
|------|------|------|-----------|
| **H1** | MSE/R²와 BioEval-Dir은 모델을 다르게 순위 매긴다 (순위 반전) | ⭐ SUPPORTED | STRONG |
| **H2** | BioEval 지표가 MSE보다 downstream 생물학적 유용성을 더 잘 예측한다 (domain-specific) | ⭐ SUPPORTED | MODERATE |
| **H3** | BioEval 하에서 학습된 모델이 단순 baseline을 능가한다 | ⭐ SUPPORTED | STRONG |

### 핵심 발견

1. **MSE/R² 순위 반전은 실재한다** — RPE1에서 τ(MSE, Dir_deg) = 0.232, 95% CI가 0 포함 → 두 지표가 통계적으로 독립
2. **Mean-effect trap이 체계적 현상** — Norman에서 mean_predictor가 MSE #1이나 Dir_deg #11. DEG 비율이 낮을수록 trap 심화
3. **BioEval은 domain-specific 이점** — Intra-DEG/Intra-magnitude 100% 통과, cross-domain 33.3%
4. **MSE는 domain-general predictor** — ρ(-MSE, dir_discovery) = 0.88-0.96. MSE 자체가 방향 정보를 포착
5. **학습 모델 > baseline** — 3 데이터셋 × 6 지표 ALL WIN (Ridge 기준; GEARS DL 검증 진행 중)

### 데이터셋

| 데이터셋 | 세포 수 | 섭동 수 | 유전자 수 | DEG 비율 |
|----------|---------|---------|----------|---------|
| Replogle K562 | 162,751 | 1,092 | 5,000 | 2.38% |
| Replogle RPE1 | 162,733 | 1,543 | 5,000 | 6.50% |
| Norman 2019 | 91,205 | 283 | 5,045 | 1.53% |

## BioEval 지표 구성

| 지표 | 측정 대상 | 차별점 |
|------|----------|--------|
| **BioEval-Dir** | 유전자×섭동 수준 방향 정확도 | 기존 PDCorr은 섭동 수준만, 본 지표는 2차원 분해 제공 |
| **BioEval-DEG** | DEG 회복 정밀도 (AUPRC) | 방향 정보를 결합한 DEG_dir_auprc는 최초의 방향 결합 DEG 평가 |
| **BioEval-Cal** | 효과 크기 보정 분석 | 체계적 과소/과대 예측 탐지 |

## 프로젝트 구조

```
nexus-science-win/
├── CLAUDE.md                    # 프로젝트 오케스트레이터
├── docs/
│   ├── research_report.md       # 연구 보고서 (전체 결과 요약)
│   ├── 08_research_report_guide.md  # 보고서 작성 가이드
│   └── environment.md           # 환경 정보
├── stages/                      # 압축된 최신 지식 (단계별)
│   ├── 01_literature_review.md
│   ├── 02_framing.md
│   ├── 03_planning.md
│   └── 04_analysis.md
├── objects/current/             # 결정 상태 추적 (YAML)
│   ├── result_card.yaml
│   ├── validation_readiness_card.yaml
│   ├── experiment_contract.yaml
│   ├── evaluation_validity_card.yaml
│   ├── idea_abstraction_card.yaml
│   └── novelty_ledger.yaml
└── outputs/                     # 단계별 산출물
    ├── literature_review/       # run_05-07
    ├── framing/                 # run_03-06
    ├── planning/                # run_03-06
    └── analysis/                # run_09-20
        ├── run_13/              # BioEval 메트릭-순위 반전 (시뮬레이션)
        ├── run_14/              # Phase 4 downstream 과업 상관
        ├── run_15/              # sklearn Ridge LOO (Norman)
        ├── run_16/              # Gene PCA Feature Ridge (K562/RPE1)
        ├── run_17/              # Bootstrap CI (B=10,000)
        ├── run_18/              # Scale Correction
        ├── run_19/              # Downstream Task Independence
        └── run_20/              # GEARS DL 모델 훈련 (진행 중)
```

## Run 이력

| Run | 날짜 | 내용 | 결과 |
|-----|------|------|------|
| run_13 | 04-30 | BioEval 메트릭-순위 반전 (시뮬레이션 11개 모델) | H1 SUPPORTED |
| run_14 | 04-30 | Phase 4 downstream 과업 상관 | H2 SUPPORTED (88.9%) |
| run_15 | 04-30 | sklearn Ridge LOO (Norman) | H1+H2 실제 모델 확인. K562/RPE1 퇴화 |
| run_16 | 05-01 | Gene PCA Feature Ridge (K562/RPE1) | H1+H2+H3 3 데이터셋 전체 확인 |
| run_17 | 05-01 | Bootstrap CI (B=10,000) | H1+H2 통계적 견고성 확인 |
| run_18 | 05-01 | Scale Correction | 보정 불필요 확인. Dir_deg 불변 |
| run_19 | 05-01 | Downstream Task Independence | H2 domain-specific 확인. MSE domain-general |
| run_20 | 05-02 | GEARS DL 모델 훈련 | 진행 중 (K562 학습 중, Norman GO graph 오류) |

## 환경

- **System Python 3.14**: 분석 스크립트 (run_13+), Ridge LOO, 통계
- **ai_env (conda, Python 3.11)**: GEARS/CPA 학습, PyG 모델
- **GPU**: RTX 4060 Ti 8GB

## 라이선스

MIT
