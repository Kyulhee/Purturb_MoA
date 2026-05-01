# Stage 01 — Literature Review

## Workflow
1. `docs/01_literature_review.md` 가이드 확인
2. 문헌 검색 및 분석 수행
3. 산출물 → `outputs/literature_review/run_XX/`에 저장
4. 아래 지식 업데이트 (검증된 인사이트만 통합)
5. `objects/current/idea_abstraction_card.yaml` 업데이트 (아이디어 추상화, 동등 표현, 검색 커버리지, 위험 플래그)
6. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 문헌 5편 이상 리뷰 완료
- SOTA 베이스라인 수치 최소 1개 확보
- 사용 가능한 데이터셋 확인 완료
- 기존 방법론의 한계(gap) 명시 완료

---

## 1. 연구 배경

Perturb-seq 섭동 예측은 20+ 신규 방법(2024-2026)에도 불구하고 평가 지표의 근본적 한계에 직면: MSE/R²/Pearson은 재구성 정확도만 측정하며, 생물학적 충실도(DEG 회복, 방향 정확도, downstream 과업 성과)를 측정하지 못한다. 이것이 Ahlmann-Eltze(2025) Nature Methods의 "DL ≤ linear baseline" 위기의 근본 원인일 가능성이 높다.

핵심 질문: **섭동 예측 평가 지표가 생물학적 충실도를 측정하지 못한다면, 생물학적 충실도를 측정하는 지표를 설계하고, 이 지표 하에서 베이스라인 위기가 해소되는가?**

## 2. SOTA

### 평가 지표 현황
| 지표 | 측정 대상 | 한계 | 비고 |
|------|----------|------|------|
| MSE/R² | 분포 매칭 | 생물학적 방향·크기 무시; mean-effect trap | 표준 |
| Pearson | 선형 상관 | 방향 구분 불가; 아웃라이어에 민감 | 표준 |
| DEG overlap | DEG 회복 | 이진 임계값으로 정보 손실; 방향 미포함 | CPA 평가 |
| PDS | 분포 거리 | 거리 지표 선택에 민감 | ARC Virtual Cell Challenge |
| AUPRC (Zhu 2025) | DEG precision-recall | DEG 식별만, 방향/보정/downstream 상관 없음 | Briefings Bioinformatics |
| PDCorr (SCALE 2026) | 섭동 방향 상관 | 유전자 수준 분해 없음, 구현 민감성 | Cell-Eval 프레임워크 |
| Shesha stability (Raju 2026) | 기하학적 안정성 | 단일 지표, 통합 프레임워크 없음 | magnitude와 0.75-0.97 상관 |

### 베이스라인 위기 증거
| 출처 | 핵심 발견 | 비고 |
|------|----------|------|
| Ahlmann-Eltze (2025) Nature Methods | DL ≤ linear baseline across 7+ benchmarks | 원인 분석 없음 |
| Csendes (2025) BMC Genomics | FM ≤ mean predictor | 낮은 섭동 특이 분산 확인 |
| SCALE (2026) | MSE가 mean-effect trap 유발 | PDCorr+DE overlap 대안 제안 |
| 우리 run_12 | CPA > FCR (0.430 vs 0.367) | MSE 기준 — BioEval에서 반전 가능? |
| 우리 run_12 | prod rho=0.437 PASS vs A7 rho=0.326 PARTIAL | 지표 선택이 결론 변경 |

### 직접/간접 경쟁자 (5개 부분)
| 경쟁자 | 핵심 기여 | 우리와의 차이 |
|--------|----------|--------------|
| Wei et al. (2026) Nature Methods | 27방법×6지표 벤치마크 | 기존 지표 사용, 새 지표 설계 없음, 순위 반전 분석 없음 |
| Zhu et al. (2025) Briefings Bioinf | AUPRC 지표 | DEG 식별만, 방향/보정/downstream 상관 없음 |
| SCALE/Chen et al. (2026) | Cell-Eval (PDCorr+DE overlap) | 구현 민감성 문제, 이론적 분석 없음 |
| Shesha/Raju (2026) | 기하학적 안정성 | 단일 지표, 통합 프레임워크 없음 |
| Csendes et al. (2025) BMC Genomics | FM ≤ mean 벤치마크 | 지표 분석 없음 |

## 3. 이론적 기반

### Goodhart's law in perturbation prediction
- 최적화 목표(평가 지표)가 모델 행동을 결정하므로, 잘못된 지표는 잘못된 모델 비교를 낳음
- MSE 최적화 → mean prediction (안전한 평균) → 생물학적 이질성 smoothing
- 생물학적 충실도를 측정하는 지표가 없으면 모델 비교가 무의미

### 3개 독립 그룹의 수렴 증거
1. Ahlmann-Eltze: DL ≤ baseline (위기 진단)
2. SCALE: MSE가 mean-effect trap 유발 (원인 지적)
3. Shesha: magnitude ≠ stability (지표 분리 필요)

## 4. Gap (Framing으로 전달)

**Gap 2: Evaluation Metrics Don't Capture What Matters ⭐ USER SELECTED**

1. **지표-생물학 상관 부재**: 기존 지표 중 어떤 것이 downstream 생물학적 유용성을 예측하는지 정량 분석 없음
2. **생물학적 충실도 지표 설계 공백**: 유전자 수준 분해능 + 방향 인식 + 보정 분석 + 구현 견고성을 모두 갖춘 지표 부재
3. **지표-순위 반전 분석 미수행**: 지표 교체가 모델 순위를 바꾸는지 누구도 정량 분석하지 않음 (핵심 차별화)
4. **베이스라인 위기 원인 불명**: 지표 아티팩트인지 실재하는지 판별 불가

---

## Run 이력 (세부 내용은 outputs/literature_review/run_XX/ 참조)
- **run_07** (2026-04-30): BioEval 신규 프로젝트 문헌 리뷰. 5개 Gap 식별, Gap 2(평가 지표)를 사용자 선택. Ahlmann-Eltze/SCALE/Shesha/Zhu/Csendes 핵심 논문 리뷰. 5개 부분 경쟁자 확인
- run_06 (2026-04-30): [FCR-ICM] Direction A 심층 리뷰
- run_05 (2026-04-29): [FCR-ICM] 5개 방향 탐색
- run_01-04 (2026-04-25~29): [FCR-ICM] 초기 문헌 리뷰
