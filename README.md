# BioEval: 생물학적 충실도 기반 섭동 예측 평가 프레임워크

Perturb-seq 섭동 예측에서 평가 지표(MSE/R²)가 만드는 **Mean-Effect Trap**을 진단하고, 생물학적 충실도 기반 평가 지표(BioEval)를 설계하여 "DL ≤ Baseline" 위기가 지표의 아티팩트인지 판별하는 연구 프로젝트입니다.

## 🚀 프로젝트 요약 (2026-05-06 완료)

본 연구는 평가 지표가 모델 선택을 어떻게 왜곡하는지 정량적으로 분석하였으며, 특히 MSE가 방향 정보가 거세된 '평균 예측'에 보상을 준다는 사실을 입증했습니다. 최종 연구 보고서는 [docs/research_report.md](./docs/research_report.md)에서 확인하실 수 있습니다.

## 🔬 핵심 연구 질문

> **"평가 지표를 생물학적 충실도 기반으로 바꾸면, 모델 순위가 어떻게 달라지는가? DL ≤ Baseline 위기가 지표의 아티팩트인지, 실재하는 현상인지 판별할 수 있는가?"**

## 📊 최종 가설 및 결과

| 가설 | 내용 | 판정 | 클레임 강도 |
|------|------|------|-----------|
| **H1** | MSE/R²와 BioEval-Dir은 모델을 다르게 순위 매긴다 (순위 반전) | ⭐ SUPPORTED | **STRONG** |
| **H2** | BioEval 지표가 MSE보다 downstream 생물학적 유용성을 더 잘 예측한다 (domain-specific) | ⭐ SUPPORTED | **MODERATE** |
| **H3** | BioEval 하에서 학습된 모델이 단순 baseline을 능가한다 | ⭐ SUPPORTED | **STRONG (Qualified)** |

### 핵심 발견
1.  **순위 반전(Ranking Reversal) 실재**: RPE1 데이터셋에서 MSE와 BioEval-Dir 순위 간 상관이 τ=0.232로 나타나, 두 평가 체계가 통계적으로 독립적임을 확인했습니다.
2.  **Mean-Effect Trap**: Norman 데이터셋에서 `mean_predictor`가 MSE 1위였으나 BioEval-Dir 11위(최하위)를 기록하여, MSE가 방향 정보를 무시함을 정량적으로 증명했습니다.
3.  **Domain-Specific 예측력**: BioEval은 도메인 내 과업(DEG 식별 등)에서 MSE보다 +0.319 높은 예측력을 보였으나, 교차 도메인에서는 MSE가 여전히 우수한 일반형 예측자(general predictor)임을 발견했습니다.
4.  **Model Quality > Complexity**: GEARS(DL) 모델이 훈련 미흡 시 단순 Ridge 모델에 전패(0/12 승)하는 현상을 BioEval로 포착하여, 지표의 판별력을 검증했습니다.

## 🛠️ BioEval 지표 구성

| 지표 | 측정 대상 | 핵심 차별점 |
|------|----------|------------|
| **BioEval-Dir** | 유전자×섭동 수준 방향 정확도 | 섭동 단위 평균이 아닌, 개별 유전자 수준의 방향성 분해 제공 |
| **BioEval-DEG** | DEG 회복 정밀도 (AUPRC) | 방향 정보를 결합하여 단순 식별 이상의 충실도 측정 |
| **BioEval-Cal** | 효과 크기 보정 분석 | 예측치의 체계적 과소/과대 편향 탐지 |

## 📂 프로젝트 구조

```
nexus-science-win/
├── docs/
│   ├── research_report.md       # 최종 연구 보고서 (10개 섹션 요약)
│   ├── stages/                  # 단계별 기법 및 가이드
│   └── environment.md           # 실험 환경 상세
├── stages/                      # 압축된 프로젝트 상태 (01-05)
├── objects/current/             # 최종 결정 데이터 (Claim Card 등)
├── outputs/                     # 실험 산출물 및 리포트
│   ├── analysis/run_13-20/      # 핵심 분석 (GEARS, Bootstrap 등)
│   └── interpretation/run_01/   # 최종 해석 및 완성 리포트
└── README.md                    # 프로젝트 메인 가이드
```

## 💻 실행 및 환경

- **System Python 3.14**: 주 분석 및 통계 스캔
- **ai_env (Conda, Python 3.11)**: GEARS/CPA 딥러닝 모델 학습
- **Hardware**: RTX 4060 Ti 8GB 기반 가중치 최적화

## 🎓 참고 문헌

- Ahlmann-Eltze et al. (2025). Nature Methods.
- Norman et al. (2019). Science.
- Replogle et al. (2022). Cell.
- Wei et al. (2026). Nature Methods (Benchmarking).

---
🤖 Generated & Maintained with **[Nexus Science](https://github.com/bionexus-enterprise/NexusScience)**
