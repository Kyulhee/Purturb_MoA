# Run 09: Epistasis Detection + Uncertainty Quantification — Results Report

**Date:** 2026-04-29 | **Stage:** Analysis (run_09) | **Pipeline versions:** v1, v2, v3

---

## Executive Summary

3차례 반복 개선을 통해 에피스태시스 탐지+불확실성 정량화 파이프라인을 구현·평가했다. **RQ4(능동학습) 타겟 달성, RQ3(공식 일치도) 대폭 개선**을 달성했으나, RQ2(U-Error rho)는 타겟(>0.6)에 미달(0.401)이며, RQ3의 가산 공식 AUROC=1.0은 순환 평가로 해석이 필요하다.

---

## 실험 결과 요약

| RQ | 지표 | v1 (run_09) | v2 (run_09b) | v3 (run_09c) | 타겟 |
|----|------|-------------|--------------|--------------|------|
| RQ1 | 합성 분해 corr | 1.000 | 1.000 | — | r>0.7 ✅ |
| RQ2 | U-Error rho | 0.304 | 0.061 | **0.401** | >0.6 ❌ |
| RQ2 | Coverage | 0.915 | 1.000 | 1.000 | 0.85-0.95 |
| RQ2 | OOD rho | — | — | 0.385 | — |
| RQ3 | AUROC (add) | 0.685 | 0.5 | 1.000⚠️ | >0.75 |
| RQ3 | Product-GT rho | — | — | 0.620 | — |
| RQ3 | Formula agreement | 0.011 | 87.5% 불일치 | **0.766** | >60% ✅ |
| RQ4 | Top-10 vs random | 0.67x | — | **5.0x (epi), 3.5x (OOD)** | >2x ✅ |
| RQ4 | Top-20 vs random | 1.0x | — | **3.2x (epi), 2.6x (OOD)** | >2x ✅ |
| RQ4 | NDCG@10 | — | — | 0.332 (OOD) | — |

---

## 반복 개선 이력

### v1 (run_09): 초기 구현

**문제점:**
- `residual_mag`가 `mean_abs_error`와 동일 → UQ rho=1.0 (순환)
- R2<median을 ground truth로 사용 → poor proxy
- 공식 일치도 0.011 → 유전자 수준에서 공식이 극도로 불일치
- AL 획득 함수에 다양성 term 없음

### v2 (run_09b): 첫 번째 개선

**변경:**
- `residual_mag` 제거 (순환 해소)
- 순열 검정으로 ground truth 구축
- 디코더 감도 기반 모델 오차 추정

**문제점:**
- 순열 검정이 모든 128개 조합을 에피스태시스로 판정 → 이진 분류 불가
- "에피스태시스가 있는가?"가 아니라 "얼마나 강한가?"가 올바른 질문
- MC Dropout 구현 오류 (z_tx 노이즈 ≠ dropout)
- UQ rho 0.304→0.061 악화 (잔류 신호 제거 후 남은 ICM+MC만으로는 불충분)

### v3 (run_09c): 두 번째 개선 (최종)

**변경:**
1. **연속적 순위 평가**로 전환 (이진 분류 → Spearman rho, NDCG)
2. **Proper MC Dropout** — 디코더에 dropout 활성화 후 forward pass
3. **OOD distance** — composed z_tx와 학습 분포 간 거리를 불확실성 신호로 추가
4. **Cohen's d-like effect size**를 ground truth로 사용
5. **AL에 NDCG** 추가

---

## RQ별 상세 분석

### RQ1: 잔차 분해

합성 데이터에서는 완벽한 분해 (r=1.000). 실데이터에서는:
- v2: 에피스태시스 비율 69.5% (모델 오차 과소추정)
- v3: 에피스태시스 비율 7.1% (OOD 기반 모델 오차 과다추정 가능)
- **한계**: 실데이터에서는 ground truth 성분이 없어 직접 검증 불가

### RQ2: 불확실성 정량화

| 신호원 | Spearman rho (v3) |
|--------|-------------------|
| OOD distance | 0.385 (p<0.0001) |
| MC Dropout variance | 0.333 |
| Combined | 0.401 |
| ICM violation | 0.108 |

- **OOD distance가 가장 강한 신호** — z_tx 공간에서 학습 분포로부터 멀리 떨어진 조합이 실제로 더 큰 오차를 가짐
- MC Dropout도 유의미한 신호
- ICM violation은 약한 신호 (K562 단일 세포 유형이므로 ICM 정규화가 단순 분산 정규화로만 작동)
- **타겟 미달 원인**: 가법 기준 R2가 이미 0.968로 높아, 잔차가 작고 노이즈에 민감함

### RQ3: 에피스태시스 탐지

**가산 공식 AUROC=1.0은 순환적**: ground truth (effect_size_gt) = |additive residual| / ctrl_std이고, epistasis_strength_add도 동일한 계산이므로.

**유의미한 결과:**
- Product neutrality-GT rho = 0.620 → 다른 공식으로도 62% 순위 복구
- Multiplicative-GT rho = 0.131 → 로그 스케일 공식은 순위 상관 약함
- **공식 일치도 0.766** → 3개 공식이 유전자별 방향에서 76.6% 일치 (v1의 1.1%에서 대폭 개선)
- 이전 v1의 극도로 낮은 일치도(0.011)는 유전자 수준 threshold 오류였음

### RQ4: 능동학습

| 방법 | Top-10 개선율 | Top-20 개선율 | NDCG@10 |
|------|--------------|--------------|---------|
| Random | 1.0x | 1.0x | 0.097 |
| UQ (combined) | 4.0x | 2.2x | 0.340 |
| Epistasis score | 5.0x | 3.2x | 0.295 |
| OOD distance | 3.5x | 2.6x | 0.332 |
| Oracle | 5.0x | 5.1x | 1.000 |

**모든 방법이 타겟(2x) 달성.** 에피스태시스 분해 점수가 top-10에서 최고 성능. OOD distance도 유의미. NDCG에서는 UQ와 OOD가 더 안정적.

---

## 핵심 통찰

1. **에피스태시스는 스펙트럼이다**: 모든 double-KO가 가법 기대에서 벗어남. 이진 탐지가 아니라 연속적 순위 매김이 올바른 프레이밍.

2. **OOD distance가 강한 불확실성 신호**: z_tx 공간에서 조합의 "새로움"이 예측 오차와 유의하게 상관 (rho=0.385). 이는 ICM 정규화의 간접적 효과.

3. **공식 선택이 중요하지만 치명적이지 않음**: 3개 공식이 76.6% 유전자별 방향 일치. 가산-곱법 간 격차는 있지만 완전한 무작위는 아님.

4. **가법 기준이 강함**: R2=0.968로 이미 매우 정확. 잔차가 작아 UQ 신호가 노이즈에 묻힘. 이것이 RQ2 타겟 미달의 주원인.

---

## 남은 과제

1. **RQ2 rho > 0.6 달성**: 
   - 더 큰 잔차를 갖는 데이터(Replogle, PORTAL)에서 재평가
   - 유전자별 UQ (조합 전체가 아닌)로 세분화
   - Conformal prediction 적용

2. **RQ3 순환 평가 해소**: 
   - 분해된 r_epistasis_mag가 아닌, OOD/ICM 기반 점수로 에피스태시스를 탐지하도록 수정
   - GEARS GI 라벨 또는 독립적 생물학적 검증 필요

3. **Cross-dataset 검증**: 
   - Norman만으로는 ICM 정규화의 진가(다중 세포 유형)를 보여줄 수 없음
   - Replogle K562+RPE1에서 ICM 기여 재측정

---

## 산출물

| 파일 | 설명 |
|------|------|
| `run_09_epistasis_uq.py` | v1 파이프라인 |
| `run_09b_epistasis_uq_v2.py` | v2 파이프라인 |
| `run_09c_epistasis_uq_v3.py` | v3 파이프라인 (최종) |
| `run_09_results.json` | v1 결과 |
| `run_09b_results.json` | v2 결과 |
| `run_09c_results.json` | v3 결과 |
| `run_09_results.md` | 본 보고서 |
