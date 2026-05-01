# BioEval Experiment Plan — Run 06

**Date:** 2026-04-30 | **Direction:** Evaluation Metrics for Perturbation Prediction (Gap 2)
**Predecessor:** Framing run_06

---

## 1. 연구 질문 (from stages/02)

**"섭동 예측 평가 지표가 생물학적 충실도를 측정하는가? 생물학적 충실도를 측정하는 지표를 설계하고, 이 지표 하에서 Ahlmann-Eltze의 'DL ≤ baseline' 위기가 해소되는가?"**

### 하위 질문
- **RQ1**: 기존 지표 중 어떤 것이 생물학적 유용성과 상관하는가?
- **RQ2**: BioEval 지표를 설계할 수 있는가? 기존 지표보다 downstream 생물학적 유용성을 더 잘 예측하는가?
- **RQ3**: BioEval 하에서 DL > baseline 체계가 존재하는가? 베이스라인 위기가 지표 아티팩트인가?

---

## 2. 방법론 아키텍처 (4-Phase)

```
Phase 1: 데이터 확보 + 모델 예측 수집
Phase 2: BioEval 지표 구현 + 기존 지표 계산
Phase 3: 지표-순위 반전 분석 (RQ1 + RQ3)
Phase 4: 지표-downstream 과업 상관 분석 (RQ2)
```

### Phase 1: 데이터 확보 + 모델 예측 수집

**목표:** 다수 모델의 예측 결과를 동일 데이터셋에서 확보

**1.1 Ahlmann-Eltze 벤치마크 예측 확보**
- Ahlmann-Eltze et al. (2025) Nature Methods 부록/Supplementary에서 모델 예측치 확인
- GitHub 코드: 모델 재현 또는 사전 계산 결과 다운로드
- 대안: 재현 불가 시, 직접 모델 학습 (아래 1.2)

**1.2 직접 모델 학습 (Ahlmann-Eltze 미확보 시)**
| 모델 | 구현 | 학습 데이터 | 비고 |
|------|------|----------|------|
| Mean predictor | `np.mean(X_ctrl, axis=0)` | 제어 세포 | 베이스라인 |
| Additive linear | `sklearn.LinearRegression` | Replogle K562 | Ahlmann-Eltze 최우수 |
| CPA | `cpa` 패키지 | Replogle K562+RPE1 | 조합 오토인코더 |
| GEARS | `gears` 패키지 | Replogle K562 | GNN+GRN |
| scGPT + linear | `scgpt` 임베딩 + linear | Replogle | 사전학습 임베딩 |

**1.3 데이터 전처리**
- Replogle 2022: K562 + RPE1, 공유 848 섭동 필터링
- Norman 2019: 128 double-KO, GI ground-truth (Norman et al. 제공)
- 정규화: log1p → 라이브러리 크기 보정 (scanpy 표준)
- DEG 정의: |logFC| > 0.25 & adjusted p < 0.05 (섭동 vs 제어)

**1.4 데이터 분할**
- Replogle: train (80%) / test (20%) — 섭동 기준 분할 (세포 기준이 아님)
- Norman: train = single KO, test = double KO (조합 일반화 평가)
- 교세포: K562 train → RPE1 test (cross-cell-type 전이)

### Phase 2: BioEval 지표 구현 + 기존 지표 계산

**2.1 기존 지표 계산 (베이스라인 비교군)**

| 지표 | 구현 | 비고 |
|------|------|------|
| MSE | `sklearn.mean_squared_error` | 플랫 벡터 비교 |
| R² | `scipy.stats.pearsonr` → r² | 유전자별 → 섭동별 집계 |
| Pearson r | `scipy.stats.pearsonr` | 유전자별 → 섭동별 집계 |
| DE overlap | Jaccard(DEG_pred, DEG_obs) | 임계값 스윕: |logFC| > {0.1, 0.25, 0.5} |
| AUPRC | `sklearn.average_precision_score` | DEG 이진 분류로 변환 |
| PDCorr (SCALE) | 섭동 방향 Pearson 상관 | 구현: SCALE Cell-Eval 참조 |
| Shesha stability | 방향 일관성 (cosine sim) | 구현: Shesha 논문 수식 참조 |

**2.2 BioEval-Dir: 유전자×섭동 수준 방향 정확도**

```python
def bioeval_dir(y_true, y_pred, gene_names, perturbations):
    """
    유전자×섭동 수준 방향 정확도 + 크기 비율
    
    y_true: (n_perturbations, n_genes) 관측 logFC
    y_pred: (n_perturbations, n_genes) 예측 logFC
    
    Returns:
        dir_accuracy: 부호 일치 비율 (0-1)
        mag_ratio: |ŷ/y| 중앙값 (1.0=완벽)
        dir_per_gene: 유전자별 방향 정확도
        dir_per_pert: 섭동별 방향 정확도
    """
    # 방향 정확도
    sign_match = (np.sign(y_pred) == np.sign(y_true)).astype(float)
    dir_accuracy = sign_match.mean()
    dir_per_gene = sign_match.mean(axis=0)  # 유전자별
    dir_per_pert = sign_match.mean(axis=1)  # 섭동별
    
    # 크기 비율 (0으로 나누기 방지)
    eps = 1e-8
    mag_ratio = np.median(np.abs(y_pred) / (np.abs(y_true) + eps))
    
    return dir_accuracy, mag_ratio, dir_per_gene, dir_per_pert
```

**임계값 민감도 테스트:**
- DEG만 평가: |y_true| > threshold ({0.1, 0.25, 0.5, 1.0})
- 모든 유전자 평가 (가중치: |y_true| 비례)
- 결과: 임계값에 따른 dir_accuracy 변화 보고

**2.3 BioEval-Cal: 효과크기 보정 분석**

```python
def bioeval_cal(y_true, y_pred, gene_names, perturbations):
    """
    logFC 예측 vs 실제의 보정 분석
    
    Returns:
        cal_slope: 회귀 기울기 (이상=1.0)
        cal_intercept: 회귀 절편
        underpredict_frac: 기울기 < 0.8인 유전자 비율
        overpredict_frac: 기울기 > 1.2인 유전자 비율
        cal_per_gene: 유전자별 회귀 기울기
    """
    # 전체 보정 곡선
    slope, intercept, r, p, se = scipy.stats.linregress(y_true.flatten(), y_pred.flatten())
    
    # 유전자별 보정
    cal_per_gene = []
    for g in range(y_true.shape[1]):
        sg, _, _, _, _ = scipy.stats.linregress(y_true[:, g], y_pred[:, g])
        cal_per_gene.append(sg)
    
    underpredict_frac = np.mean(np.array(cal_per_gene) < 0.8)
    overpredict_frac = np.mean(np.array(cal_per_gene) > 1.2)
    
    return slope, intercept, underpredict_frac, overpredict_frac, cal_per_gene
```

**2.4 BioEval-DEG: DEG precision-recall 곡선**

```python
def bioeval_deg(y_true, y_pred, thresholds=np.arange(0.1, 2.0, 0.1)):
    """
    DEG precision-recall 곡선 (AUPRC + 곡선 전체)
    방향 결합: sign mismatch 시 penalty
    
    Returns:
        auprc: 곡선 아래 면적
        pr_curve: (precision, recall) 쌍 리스트
        dir_aware_auprc: 방향 인식 AUPRC
    """
    # 표준 AUPRC
    y_true_bin = (np.abs(y_true) > 0.25).astype(int)
    y_pred_score = np.abs(y_pred)
    auprc = average_precision_score(y_true_bin.flatten(), y_pred_score.flatten())
    
    # 방향 인식: 부호가 틀린 DEG는 penalty
    sign_match = (np.sign(y_pred) == np.sign(y_true))
    dir_aware_score = y_pred_score * sign_match  # 방향 틀리면 0
    dir_aware_auprc = average_precision_score(y_true_bin.flatten(), dir_aware_score.flatten())
    
    return auprc, dir_aware_auprc
```

**2.5 BioEval-Composite**

```python
def bioeval_composite(dir_acc, cal_slope, deg_auprc, 
                       w_dir=0.4, w_cal=0.3, w_deg=0.3):
    """
    통합 지표: Dir + Cal + DEG 가중합
    가중치는 RQ1 결과에서 학습 (downstream 과업 상관 기준)
    """
    cal_score = 1.0 - abs(cal_slope - 1.0)  # 기울기 1.0에서 최대
    return w_dir * dir_acc + w_cal * cal_score + w_deg * deg_auprc
```

### Phase 3: 지표-순위 반전 분석 (RQ1 + RQ3)

**3.1 모델 순위 산출**

각 지표(MSE, R², Pearson, DE overlap, AUPRC, PDCorr, Shesha, BioEval-Dir, BioEval-Cal, BioEval-Composite)에 대해:
1. 각 모델의 지표 점수 계산 (섭동별 → 평균)
2. 지표 점수 기준 모델 순위 산출
3. 순위 간 비교

**3.2 Kendall τ 순위 상관 분석**

```python
from scipy.stats import kendalltau

# MSE 순위 vs BioEval 순위
tau_mse_dir, p_mse_dir = kendalltau(rank_mse, rank_bioeval_dir)
tau_mse_cal, p_mse_cal = kendalltau(rank_mse, rank_bioeval_cal)
tau_mse_comp, p_mse_comp = kendalltau(rank_mse, rank_bioeval_composite)

# 해석:
# τ > 0.7: 순위 유지 (지표 변경이 모델 선택에 영향 없음)
# 0.5 ≤ τ ≤ 0.7: 부분 반전
# τ < 0.5: 순위 반전 (지표가 모델 선택을 바꿈)
```

**3.3 순위 변동 체계적 분석**

| 분석 차원 | 질문 | 방법 |
|-----------|------|------|
| 섭동 유형별 | synergistic GI에서만 DL이 우위? | Norman GI subtype별 τ 분석 |
| 세포유형별 | 특정 세포유형에서만 반전? | K562 vs RPE1 각각 τ 계산 |
| 유전자 특성별 | 고발현/저발현 유전자에서 차이? | 유전자 발현량 사분위별 dir_accuracy |
| 효과크기별 | 큰 효과 vs 작은 효과에서 차이? | |logFC| 사분위별 분석 |

**3.4 베이스라인 위기 원인 판별 (RQ3 핵심)**

| 시나리오 | Kendall τ | 해석 | 임팩트 |
|----------|-----------|------|--------|
| 순위 반전 | < 0.5 | BioEval이 MSE와 다른 모델을 선택 → 위기는 지표 아티팩트 | 높음 |
| 순위 유지 | > 0.7 | BioEval 하에서도 baseline이 우위 → 위기 실재 | 중간 |
| 부분 반전 | 0.5-0.7 | 특정 체계에서만 DL 우위 → 미묘한 결과 | 높음 |

### Phase 4: 지표-downstream 과업 상관 분석 (RQ2)

**4.1 downstream 과업 정의**

| 과업 | 정의 | 평가 방법 |
|------|------|----------|
| DEG 회복 | top-k DEG 식별 정확도 | Precision@k, Recall@k (k=10,20,50) |
| AL 효율 | 지표로 선택된 모델이 AL에서 더 나은가? | 모의 AL: 불확실성 기반 순차 선택 vs 랜덤 |
| Hit prioritization | 지표 순위가 실제 hit 순위와 일치하는가? | Spearman ρ (지표 순위 vs hit 순위) |
| 교세포 전이 | 지표로 선택된 모델이 다른 세포유형에서 더 나은가? | Cross-CT R² (K562 선택 → RPE1 성과) |

**4.2 지표-과업 상관 계산**

```python
# 각 지표의 모델 순위 vs 각 downstream 과업의 모델 순위
for metric_name in all_metrics:
    for task_name in downstream_tasks:
        rho, p = spearmanr(rank_metric[metric_name], rank_task[task_name])
        results.append({
            'metric': metric_name,
            'task': task_name,
            'spearman_rho': rho,
            'p_value': p
        })

# BioEval이 MSE보다 downstream 과업을 더 잘 예측하는지 비교
# 타겟: rho(BioEval, task) > rho(MSE, task) by ≥ 0.1
```

**4.3 BioEval-Composite 가중치 학습**

- Grid search: w_dir ∈ {0.2, 0.3, 0.4, 0.5}, w_cal ∈ {0.1, 0.2, 0.3, 0.4}, w_deg = 1 - w_dir - w_cal
- 기준: downstream 과업 상관의 평균 Spearman ρ 최대화
- 과적합 방지: 5-fold 교차 검증 (섭동 기준 분할)

---

## 3. 평가 타겟

| RQ | 지표 | Baseline | 타겟 | 근거 |
|----|------|----------|------|------|
| RQ1 | 지표-AL 상관 (Spearman) | 0 (MSE-AL 상관) | > 0.5 | MSE가 AL을 예측 못하면 BioEval이 해야 |
| RQ1 | 지표-DEG 상관 (Spearman) | AUPRC 기준치 | BioEval > MSE by ≥0.1 | Zhu: R²-AUPRC 불일치 입증 |
| RQ2 | BioEval-Dir 방향 정확도 | 0.5 (우연) | > 0.7 | 70%+ 방향 일치가 유용성 최소 기준 |
| RQ2 | BioEval-Cal 보정 기울기 | — | 0.8-1.2 범위 | 과소/과대 예측 20% 이내 |
| RQ3 | Kendall τ (MSE vs BioEval) | 1.0 (동일) | < 0.5 또는 > 0.7 | 이진 결과 — 둘 다 의미 있음 |
| RQ3 | Synergistic GI에서 DL>baseline | 0% (MSE 하) | > 30% (BioEval 하) | 특정 체계에서 DL 우위 탐지 |

---

## 4. 소거/민감도 실험

| ID | 변인 | 목적 |
|----|------|------|
| S1 | DEG 임계값 스윕 (0.1, 0.25, 0.5, 1.0) | BioEval-Dir 임계값 민감도 |
| S2 | BioEval-Composite 가중치 변동 | 통합 지표 견고성 |
| S3 | 방향 평가 대상: 전체 vs DEG-only | 방향 정확도의 유전자 선택 효과 |
| S4 | 보정 분석: 유전자별 vs 섭동별 vs 전체 | 보정의 분해 수준 효과 |
| S5 | 데이터셋 교차: Replogle↔Norman | 지표 일반성 |
| S6 | Baseline 모델 추가/제거 | 순위 안정성 |
| S7 | AUPRC(Zhu) vs BioEval-DEG 비교 | 기존 DEG 지표와의 직접 비교 |

---

## 5. 데이터

| 데이터셋 | 용도 | Phase | 비고 |
|----------|------|-------|------|
| Replogle 2022 | RQ1-3 주 평가 | 1-4 | K562+RPE1, 848 공유 섭동 |
| Norman 2019 | RQ1-3 조합 평가 | 1-4 | 128 double-KO, GI ground-truth |
| Ahlmann-Eltze 벤치마크 | RQ3 직접 재현 | 1 | 7+ 벤치마크, 예측 확보 시도 |
| PBMC (Zhu 2025) | RQ1 AUPRC 비교 | 3 | 7 cell types, IFN-γ |
| PORTAL 2026 | RQ3 대규모 검증 | 3 (선택) | 665K pairwise |

---

## 6. 핵심 리스크와 대응

| 리스크 | 확률 | 영향 | 대응 |
|--------|------|------|------|
| Ahlmann-Eltze 예측 확보 불가 | 중간 | 높음 | 직접 모델 학습(Phase 1.2)으로 대체. Replogle+Norman으로 재현 |
| 순위 반전이 안 일어남 | 중간 | 낮음 | 부정 결과도 발표 가능. "위기 실재" = 중요 발견 |
| BioEval 구현 민감성 | 중간 | 중간 | S1-S4 민감도 실험으로 견고성 입증 |
| Wei et al.와 중복 | 낮음 | 높음 | Wei는 기존 지표 벤치마크, 우리는 새 지표 설계+순위 반전. 본질 차이 |
| DEG 임계값이 결과 지배 | 중간 | 중간 | S1 임계값 스윕; 결과가 임계값에 견고해야 의미 |
| 단일 아키텍처 한계 | 낮음 | 중간 | 다수 모델(linear, CPA, GEARS, scGPT+linear) 비교 |

---

## 7. 이전 실패에서의 설계 원칙

| 교훈 | 출처 | BioEval 설계 반영 |
|------|------|-------------------|
| 지표 선택이 결론 변경 | run_12 (prod rho=0.437 vs A7 rho=0.326) | 이것이 RQ의 직접적 동기. 이번에는 체계적 분석 |
| 경쟁자 조기 확인 | FCR-ICM 8 runs 후 BuDDI/C3TL 발견 | Framing에서 5개 부분 경쟁자 이미 확인 |
| 단순 모델 > 복잡 모델 | CPA > FCR (0.430 vs 0.367) | MSE 기준이었음 — BioEval에서 반전 가능? |
| 잠재공간-유전자공간 갭 | run_05 | 모든 지표를 유전자 공간에서 계산 |
| 소거실험 필수 | run_09 | S1-S7 민감도/소거 매트릭스 설계 |
| AUROC=1.0 동어반복 | run_09 | Kendall τ(순위 상관)로 평가 — 순환 불가 |

---

## 8. 구현 순서

| Step | 작업 | 산출물 | 의존성 |
|------|------|--------|--------|
| 1 | Ahlmann-Eltze 코드/예측 확보 시도 | 데이터 접근 가능 여부 판정 | 없음 |
| 2 | Replogle + Norman 데이터 로드/전처리 | 전처리된 AnnData | 없음 |
| 3 | 모델 학습 (Mean, Linear, CPA, GEARS) | 모델 예측치 행렬 | Step 2 |
| 4 | 기존 지표 계산 (MSE, R², Pearson, DE overlap, AUPRC, PDCorr) | 지표 점수 테이블 | Step 3 |
| 5 | BioEval 지표 구현 (Dir, Cal, DEG, Composite) | BioEval 점수 테이블 | Step 3 |
| 6 | 민감도 실험 S1-S4 | 임계값/가중치 민감도 결과 | Step 5 |
| 7 | 모델 순위 산출 + Kendall τ 분석 | 순위 비교 테이블, τ 값 | Step 4, 5 |
| 8 | 순위 변동 체계적 분석 (섭동/유전자/효과크기별) | 분해 분석 결과 | Step 7 |
| 9 | downstream 과업 정의 + 평가 | downstream 과업 점수 테이블 | Step 3 |
| 10 | 지표-downstream 상관 분석 | Spearman ρ 행렬 | Step 4, 5, 9 |
| 11 | 교차 데이터셋 검증 (Replogle↔Norman) | 일반성 결과 | Step 7, 10 |
| 12 | 결과 종합 + 해석 | result_card.yaml | Step 6-11 |
