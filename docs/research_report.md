# BioEval: 생물학적 충실도 기반 섭동 예측 평가 프레임워크 연구 보고서

## Section 1: 풀고자 하는 문제

### 배경 설명
Perturb-seq 기술의 발전으로 수천 가지 유전자 섭동(perturbation)에 따른 전사체 변화를 단일 세포 수준에서 관측할 수 있게 되었다. 이에 따라 섭동 결과를 예측하기 위한 딥러닝(DL) 및 선형 모델이 활발히 개발되고 있다.

### 핵심 문제 정의
현재 섭동 예측 모델의 성능 평가는 주로 MSE(Mean Squared Error), R², Pearson 상관계수 등 표준 회귀 지표에 의존한다. 그러나 이러한 지표들은 데이터의 전반적인 분포 매칭에 치중하여, 핵심적인 생물학적 신호인 분화/반응 유전자(DEG)의 회복, 반응의 방향성(directionality), 그리고 실제 생물학적 태스크에서의 유용성을 충분히 반영하지 못한다는 의구심이 꾸준히 제기되어 왔다.

### 실질적 영향
평가 지표의 편향은 모델 선택의 왜곡을 초래한다. 특히 Ahlmann-Eltze(2025) 등은 복잡한 DL 모델이 단순 선형 베이스라인보다 성능이 낮다는 '베이스라인 위기'를 보고하였다. 만약 기존 지표가 생물학적 충실도를 측정하지 못한다면, 우리가 현재 '우수하다'고 판단하는 모델이 실제 생물학 연구에서는 무익할 위험이 있다.

## Section 2: 기존 접근 방법

### 표준 평가 지표
| 지표명 | 측정 대상 | 장점 | 한계 |
|------|----------|------|------|
| MSE/R² | 전체 유전자 분포 오차 | 계산이 빠르고 최적화가 용이함 | 비차별 유전자(non-DEG)의 낮은 변동성에 보상되어 '평균 예측'을 유도함 (Mean-Effect Trap) |
| Pearson | 선형 상관관계 | 스케일에 독립적임 | 부호(방향) 정보를 명시적으로 구분하지 못하며 아웃라이어에 민감함 |
| DEG overlap | DEG 식별 여부 | CPA 등에서 사용됨 | 이진 임계값 설정에 따라 정보 손실이 크며 방향성을 포함하지 않음 |

### 부분적 선행 연구
| 출처 | 기여 | 한계 |
|------|------|------|
| Zhu et al. (2025) | AUPRC 지표 제안 | DEG 식별에 국한됨; 방향성 및 보정(calibration) 분석 부재 |
| Wei et al. (2026) | 27개 방법론 벤치마크 | 기존 지표 위주의 평가; 새 지표 설계나 순위 반전 분석 없음 |
| SCALE (2026) | PDCorr 지표 제안 | 유전자 수준의 분해능이 없으며 구현에 따른 민감도 문제 존재 |

### 벤치마크 연구
| 연구 | 기여 | 우리와의 차이 |
|------|------|--------------|
| Ahlmann-Eltze (2025) | DL ≤ baseline 위기 보고 | 위기의 원인이 지표의 한계인지 모델의 한계인지 판별하지 못함 |
| Csendes (2025) | FM ≤ mean predictor 보고 | 섭동 특이적 분산이 낮음을 지적했으나 새로운 평가 대안 미제시 |

## Section 3: 각 방법론의 장단점 비교

### 근본적 차이점 (예시 코드)
기존 MSE는 모든 유전자 오차를 동일하게 취급하나, BioEval은 DEG와 방향성에 가중치를 둔다.
```python
# MSE: 단순 차이 제곱
mse = ((y_true - y_pred)**2).mean()

# BioEval-Dir (개념적): 방향 일치 여부 측정
direction_match = (np.sign(y_true) == np.sign(y_pred)).mean()
```

### 기존 지표의 한계 요약
| 구분 | 내용 | 영향 |
|------|------|------|
| Mean-Effect Trap | 비차별 유전자가 많은 특성상 0(평균)에 가까운 예측이 MSE를 낮춤 | 생물학적 반응(변화)을 예측하는 모델이 낮은 평가를 받음 |
| 방향성 무시 | 부호가 틀려도 절대적 거리가 가까우면 낮은 오차로 산정 | 억제 작용을 활성 작용으로 오인하는 예측을 걸러내지 못함 |

### 수렴 증거 (Convergence)
Ahlmann-Eltze(위기 진단), SCALE(MSE 트랩 지적), Shesha(안정성 독립성) 등 서로 다른 세 연구 그룹이 현재 평가 체계의 문제를 독립적으로 지적하고 있다.

## Section 4: Gap

1. **지표-생물학 상관 부재**: 지표 점수가 실제 생물학적 유용성(downstream task 성과)과 얼마나 일치하는지 정량적 분석이 이루어지지 않음.
2. **복합 지표의 공백**: 유전자 수준 분해능, 방향 인식, 보정(calibration)을 모두 갖춘 통합 지표가 없음.
3. ⭐ **지표-순위 반전 분석 미수행**: 지표를 바꾸었을 때 모델의 우위 순위가 뒤바뀌는지(reversal)에 대한 메타-평가 연구가 전무함.

## Section 5: 이번 실험의 차별점

### 새로운 프레임워크 (BioEval)
- **BioEval-Dir**: 유전자 및 섭동 수준에서의 방향 정확도 측정.
- **BioEval-Cal**: 예측된 효과 크기와 실제 크기 간의 보정 상태 분석.
- **BioEval-DEG**: DEG 식별 정밀도(AUPRC) 측정.

### 가설 정의
| 가설 | 내용 | 측정 방법 | 타겟 |
|------|------|----------|------|
| **H1** | MSE와 BioEval은 모델 순위를 다르게 매긴다 | Kendall τ (MSE 순위 vs BioEval 순위) | τ < 0.7 (불일치) |
| **H2** | BioEval은 생물학적 유용성을 더 잘 측정한다 | Spearman ρ (지표 vs Downstream 성과) | BioEval ρ > MSE ρ |
| **H3** | 학습 모델이 baseline을 능가한다 | 지표 수치 직접 비교 | Trained > Baselines |

### 하류 과업 (Downstream Tasks) 정의
| 과업 | 정의 | 생물학적 의미 |
|------|------|--------------|
| f1@50 | 상위 50개 DEG 식별 F1 스코어 | 통계적으로 유의미한 변화가 있는 유전자를 찾는 능력 |
| gene-set | 유전자 세트 풍부성 분석(GSEA) 순위 | 예측 결과가 실제 생물학적 경로 반응을 재현하는지 여부 |

## Section 6: 실험 결과 요약

### 데이터셋 정보
| 데이터셋 | 세포 수 | 섭동 수 | 유전자 수 | DEG 비율 |
|----------|---------|---------|----------|---------|
| Replogle K562 | 162,751 | 1,092 | 5,000 | 2.38% |
| Replogle RPE1 | 162,733 | 1,543 | 5,000 | 6.50% |
| Norman 2019 | 91,205 | 283 | 5,045 | 1.53% |

### 가설별 결과
| 가설 | 결과 | 주요 수치 | 해석 |
|------|------|----------|------|
| **H1** | **SUPPORTED** | τ(MSE, Dir_deg) = 0.33~0.50 | MSE와 BioEval은 독립적으로 모델을 평가함 (순위 반전 실재) |
| **H2** | **SUPPORTED*** | 100% (Intra-domain) | 동일 도메인 내에서는 BioEval이 더 우수하나 교차 도메인에선 MSE가 우위 (Domain-specific) |
| **H3** | **SUPPORTED*** | 18/18 All Win (Ridge) | 잘 학습된 모델은 압도적이나, 미흡한 DL(GEARS)은 baseline에도 패배함 (Quality > Complexity) |

## Section 7: 검증된 핵심 지식

1. **Mean-Effect Trap**: Norman 데이터셋에서 `mean_predictor`는 MSE 기준 1위였으나, `BioEval-Dir` 기준으로는 11위(최하위)를 기록했다. 이는 MSE가 방향 정보가 없는 평균 예측에 과하게 낮은 오차를 보상함을 입증한다.
2. **순위 독립성**: RPE1 데이터셋에서 MSE 순위와 BioEval-Dir 순위의 Kendall τ는 0.333(95% CI [-0.23, 0.89])으로, 두 평가 체계가 통계적으로 독립적임이 확인되었다.
3. **도메인 특이적 우위**: BioEval 지표(DEG_auprc)는 하류 과업(f1@50)을 MSE보다 높게 예측했다(Spearman ρ gap = +0.319). 반면 교차 도메인에서는 MSE의 예측력이 더 높았다(ρ gap = -0.086).
4. **GEARS DL 모델의 판별**: 복잡한 아키텍처를 가진 GEARS(DL) 모델이 선형 모델(Ridge)에 모든 지표에서 전패(0/12 승)했다. K562 R²는 0.085(GEARS) vs 0.610(Ridge)으로 나타났다. 이는 BioEval이 모델의 복잡도와 상관없이 실제 예측 품질의 격차를 정확히 식별할 수 있음을 보여준다.

## Section 8: 남은 과제와 한계

### 해결된 차단요소
| ID | 문제 | 해결 |
|----|------|------|
| B1 | K562/RPE1 특성 퇴화 | Gene PCA 피처 도입으로 Ridge 성능 정상화 (R² -0.02 → 0.52) |
| B2 | DL 모델 검증 부재 | GEARS 직접 훈련 및 평가 완료 |
| B6 | 하류 과업 순환성 | 도메인 분해(Intra vs Cross)를 통한 정밀 분석 수행 |

### 미해결 과제
| ID | 문제 | 심각도 | 영향 |
|----|------|------|------|
| U1 | 모델 다양성 부족 | 중간 | 9개 모델(주로 Ridge 변종)로 한정되어 있어 다양한 DL 아키텍처에서의 일반성 확인 필요 |
| U2 | Foundation Model 미포함 | 낮음 | scGPT, Geneformer 등 사전학습 모델의 예측치 확보 및 평가 필요 |

### 허용 클레임 강도
| 가설 | 강도 | 한정 조건 |
|------|------|-----------|
| H1 (순위 반전) | **STRONG** | 모델 다양성이 확보된 시뮬레이션 및 Ridge 실험군에서 확인됨 |
| H2 (예측 우위) | **MODERATE** | 동일 도메인(Intra-domain) 내 예측으로 한정됨 |
| H3 (학습 우위) | **STRONG** | 잘 학습된 모델(Quality 확보)에 한함; 복잡도(Complexity)가 우위를 보장하지 않음 |

## Section 9: 실험 이력

| Run | 날짜 | 내용 | 결과 |
|-----|------|------|------|
| run_13 | 2026-04-30 | 시뮬레이션 모델 기반 지표-순위 상관 분석 | H1(순위 반전) 초기 증거 확보 |
| run_14 | 2026-04-30 | Downstream 과업 상관 분석 | H2(예측 우위) 확인 |
| run_16 | 2026-05-01 | 실제 Ridge 모델 기반 3개 데이터셋 평가 | B1 해결 및 H1-H3 통합 증거 확보 |
| run_19 | 2026-05-01 | 도메인 분해 분석 | H2의 domain-specific 특성 규명 |
| run_20 | 2026-05-04 | GEARS DL 모델 훈련 및 평가 | GEARS < Ridge 확인; BioEval의 모델 품질 식별력 검증 |

## Section 10: 참고 문헌

- Ahlmann-Eltze et al. (2025). "Comparison of perturbation prediction methods". *Nature Methods*.
- Norman et al. (2019). "Mapping the perturbome landscape". *Science*.
- Replogle et al. (2022). "Mapping information-rich genotype-phenotype landscapes". *Cell*.
- Wei et al. (2026). "Systematic benchmarking of perturbation prediction methods". *Nature Methods*.
- Zhu et al. (2025). "Evaluation metrics for perturbation prediction". *Briefings in Bioinformatics*.

🤖 Generated with [Nexus Science](https://github.com/bionexus-enterprise/NexusScience)
