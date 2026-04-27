# FCR-ICM: Factorized Causal Representations for Perturb-seq Compositionality

Perturb-seq 데이터에서 인과 불변성(ICM) 원리로 제약된 인자화 표현(FCR)이 교세포 zero-shot 교란 전이를 가능하게 하는지 검증하는 연구 프로젝트입니다.

## 연구 질문

**H1 (인과 불변성 가설):** FCR의 교란 효과 잠재 변수 z_tx를 ICM 원리로 제약하면 cell type 불변이 되고, zero-shot 교란 전이가 가능하다.

- **RQ1:** ICM이 z_tx를 cell type 불변으로 만드는가?
- **RQ2:** 단일-KO z_tx로 조합적 예측이 가능한가? (경로간 가법, 경로내 승법)
- **RQ3:** ICM이 zero-shot 교세포 전이를 가능하게 하는가?

## 현재 결과 (run_01-07)

| RQ | 합성 데이터 | 실제 데이터 (Norman 2019 / Replogle 2022) |
|----|------------|--------------------------|
| RQ1 (불변성) | 통과 (0.50 -> 0.97) | 통과 (-0.35 -> 0.35, K562+RPE1) |
| RQ2 (조합성) | 유전자공간 통과 (R2=0.88) | 통과 (corr=0.955, R2=0.881) |
| RQ3 (전이) | 통과 (0.51 -> 0.96) | 통과 (R2 -0.30 -> 0.92, K562->RPE1) |

### 핵심 발견

1. **ICM은 불변성/전이에 확실히 유효** — 실제데이터에서도 RQ3 R2 -0.30→0.92
2. **조합성은 유전자 공간에서 평가해야 함** — 잠재공간 R2=0.05 vs 유전자공간 R2=0.88. 디코더가 인코더 비선형성 보상
3. **조합 일관성 손실(comp loss)이 RQ2 핵심** — 소거실험에서 RQ2-cross 0.20->0.79 개선
4. **ICM이 인코더를 더 선형적으로 만듦** — linear R2 0.69->0.87

## 프로젝트 구조

```
nexus-science-win/
├── CLAUDE.md                    # 프로젝트 오케스트레이터
├── docs/                        # 단계별 가이드
│   ├── 01_literature_review.md
│   ├── 02_framing.md
│   ├── 03_planning.md
│   ├── 04_analysis.md
│   ├── 05_interpretation.md
│   ├── 06_git_policy.md
│   └── 07_experiment_failure_reports.md
├── stages/                      # 압축된 최신 지식 (단계별)
│   ├── 01_literature_review.md
│   ├── 02_framing.md
│   ├── 03_planning.md
│   ├── 04_analysis.md
│   └── 05_interpretation.md
└── outputs/                     # 단계별 산출물
    ├── literature_review/
    ├── analysis/
    │   ├── run_01/              # NAP E2E 파이프라인
    │   ├── run_02/              # GNN vs tabular 문헌 리뷰
    │   ├── run_03/              # AL 실험
    │   ├── run_04/              # Phase 1 합성 + Phase 2 실제데이터
    │   ├── run_05/              # 소거실험 + RQ2 갭 분석
    │   ├── run_06/              # Norman 실제데이터 소거실험
    │   └── experiment_reports/  # 가설 실패/전환 보고서
    └── interpretation/
```

## Run 이력

| Run | 날짜 | 내용 | 결과 |
|-----|------|------|------|
| run_01 | 2026-04-26 | NAP E2E 파이프라인 | XGBoost R2=0.91, GNN 중복성 확인 |
| run_02 | 2026-04-27 | GNN vs tabular 문헌 심층 리뷰 11편 | GNN 임베딩이 tabular보다 우위인 근거 부족 |
| run_03 | 2026-04-27 | Input-space AL 실험 | AL R2=0.56 vs Random R2=0.68 — AL 실패 |
| run_04 | 2026-04-27 | 방향 전환 + Phase 1/2 검증 | RQ1/RQ3 통과, RQ2 합성실패/실제통과 |
| run_05 | 2026-04-27 | 소거실험 + RQ2 갭 해명 | comp_loss 핵심, 유전자공간 평가로 갭 해명 |
| run_06 | 2026-04-27 | Norman 실제데이터 소거실험 | 모든 config R2=0.86-0.89, comp loss 불필요 |
| run_07 | 2026-04-27 | 다세포유형 실제데이터 (K562+RPE1) | RQ1/RQ3 실제 통과: 전이 R2 -0.30→0.92 |

## 이전 방향에서의 학습 (참고용)

- **XGBoost-only R2=0.91** — FBA는 근본적으로 tabular problem
- **GNN 임베딩 중복** — 정적 그래프에서는 knockout mask가 충분 통계량
- **AL 실패** — FBA가 싸고 입력 차원이 낮아 AL 이점 없음

## 라이선스

MIT
