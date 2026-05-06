# Stage 04 — Analysis

## Workflow
1. `docs/04_analysis.md` 가이드 확인
2. `stages/03_planning.md`에서 실험 설계 및 기준 확인
3. `objects/current/experiment_contract.yaml` + `evaluation_validity_card.yaml` 사전 확인
4. 데이터 전처리, 모델 학습, 평가 수행
5. 결과를 Planning 타겟과 비교
6. 산출물 → `outputs/analysis/run_XX/`에 저장
7. 아래 지식 업데이트
8. `objects/current/result_card.yaml` 업데이트 (주요 지표, 베이스라인 비교, 소거 요약, 다음 단계)
9. `objects/current/validation_readiness_card.yaml` 업데이트 (증거 상태, 차단요소, 허용 클레임 강도)
10. 방향 전환 시: `objects/current/pivot_diagnosis_card.yaml` 작성 필수
11. CLAUDE.md의 `current_stage` 업데이트

## Done when
- Planning 타겟 성능 달성 또는 미달 시 사용자 보고 후 방향 결정

---
## Stage 04 → 05 전환 완료
모든 가설(H1 STRONG, H2 MODERATE, H3 STRONG qualified) 충분한 증거 확보. claim_card.yaml 작성 완료.

---

## 검증된 핵심 지식

### H1: 순위 반전 실재 (STRONG)

**핵심 결과**: MSE/R²와 BioEval-Dir은 모델을 독립적으로 순위 매긴다.

| 데이터셋 | τ(MSE, Dir_all) | τ(MSE, Dir_deg) | τ(Pearson, Dir_deg) | Bootstrap 95% CI |
|----------|:---------------:|:---------------:|:-------------------:|:----------------:|
| K562 | 0.500 | 0.611 | 0.056 | CI가 0 미포함 (중등도 불일치) |
| RPE1 | 0.333 | 0.389 | 0.111 | CI가 0 포함 (독립 확인) |
| Norman | 0.500 | 0.500 | -0.167 | CI가 0 포함 (독립 확인) |

- RPE1이 가장 명확한 H1 증거 (|τ|=0.232, 95% CI가 0 포함)
- Norman: mean_predictor가 MSE #1이나 Dir_deg #11 — Mean-Effect Trap 직접 증거
- DEG 비율이 낮을수록 Trap 심화 (Norman 1.53% < K562 2.38% < RPE1 6.50%)
- 교세포 일관성: K562↔RPE1 모든 지표 τ > 0.78

### H2: BioEval은 domain-specific 이점 (MODERATE)

**핵심 결과**: BioEval 지표가 MSE보다 downstream 생물학적 유용성을 더 잘 예측한다. 단, intra-domain에서만.

| 도메인 | 통과율 | 평균 gap | 해석 |
|:------:|:------:|:--------:|------|
| intra-DEG (DEG_auprc ↔ f1@50) | 9/9 (100%) | +0.319 | DEG_auprc가 f1@50을 잘 예측 |
| intra-magnitude (mag_rank ↔ gene-set) | 6/6 (100%) | +0.303 | 크기 정보가 gene-set 과업을 잘 예측 |
| intra-direction (Dir ↔ dir_discovery) | 3/6 (50%) | -0.013 | 부분 순환, 약한 이점 |
| cross-domain (Dir ↔ gene-set) | 8/24 (33.3%) | -0.086 | Dir이 gene-set을 MSE보다 잘 예측하지 못함 |

- MSE 자체가 domain-general predictor: ρ(-MSE, dir_discovery) = 0.88-0.96
- DEG_auprc가 가장 견고한 H2 지표 — Bootstrap CI에서 모든 데이터셋 유의
- 방향 독립적 과업(mag_rank, top100_overlap) 69.7% 통과 — 방향 순환에 의해 주도되지 않음

### H3: 학습 모델 > baseline (STRONG, qualified)

**핵심 결과**: Well-trained Ridge > baselines. GEARS(DL) < Ridge — model quality > model complexity.

**Ridge vs baselines**: 3 데이터셋 × 6 지표 = 18/18 ALL WIN

| 데이터셋 | Dir_deg (학습) | Dir_deg (baseline) | R² (학습) | R² (baseline) |
|----------|:-------------:|:-------------------:|:---------:|:-------------:|
| K562 | 0.985 | 0.606 | 0.523 | -0.013 |
| RPE1 | 0.989 | 0.666 | 0.652 | -0.013 |
| Norman | 0.986 | 0.571 | 0.643 | -0.002 |

**GEARS(DL) vs Ridge**: 0/12 승. K562 GEARS R²=0.085, RPE1 0.147, Norman -0.699. GEARS vs mean_predictor는 11/12 승.

### 모델 구성 (9개)

| 모델 | 유형 | 설명 |
|------|------|------|
| ridge / ridge_med / ridge_strong | 학습 | Ridge (alpha=1/10/100, analytical LOO) |
| noisy_ridge / sign_flip_ridge | 퇴화 | ridge + 15% 노이즈/부호반전 |
| mean_predictor / mean_effect | 베이스라인 | 관측치 평균 / 섭동 평균 효과 |
| constant_shrink / half_signal | 오라클 | true × 0.15 / 0.5 (H3에서 제외) |

### 데이터셋

| 데이터셋 | 세포 수 | 섭동 수 | 유전자 수 | DEG 비율 |
|----------|---------|---------|----------|---------|
| Replogle K562 | 162,751 | 1,092 | 5,000 | 2.38% |
| Replogle RPE1 | 162,733 | 1,543 | 5,000 | 6.50% |
| Norman 2019 | 91,205 | 283 | 5,045 | 1.53% |

---

## 해결된 차단요소

| ID | 문제 | 해결 |
|----|------|------|
| B1 | K562/RPE1 one-hot LOO 퇴화 | Gene PCA features (run_16). R² -0.027→0.523 |
| B2 | 실제 DL 모델 예측 없음 | GEARS 훈련 (run_20). GEARS < Ridge (0/12) |
| B3 | Norman logFC 스케일 불일치 | 보정 불필요 확인 (run_18). Dir_deg 불변 |
| B4 | Oracle baselines이 H3 왜곡 | H3 분석에서 제외 |
| B5 | Bootstrap CI 부재 | B=10,000 bootstrap (run_17) |
| B6 | downstream 과업 순환성 | 도메인 분해 (run_19). H2 domain-specific |

---

## Run 이력 (세부 내용은 outputs/analysis/run_XX/ 참조)
- **run_13** (04-30): BioEval 순위 반전 (시뮬레이션 11개 모델) → H1 SUPPORTED
- **run_14** (04-30): Phase 4 downstream 과업 상관 → H2 SUPPORTED (88.9%)
- **run_15** (04-30): sklearn Ridge LOO (Norman) → H1+H2 실제 모델 확인, K562/RPE1 퇴화
- **run_16** (05-01): Gene PCA Feature Ridge (K562/RPE1) → B1 해결, H1+H2+H3 전체 확인
- **run_17** (05-01): Bootstrap CI (B=10,000) → H1+H2 통계적 견고성 확인
- **run_18** (05-01): Scale Correction → 보정 불필요, B3 해결
- **run_19** (05-01): Downstream Task Independence → H2 domain-specific, B6 해결
- **run_20** (05-04): GEARS DL 모델 훈련+평가 → GEARS < Ridge, B2 해결
