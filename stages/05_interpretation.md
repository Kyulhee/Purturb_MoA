# Stage 05 — Interpretation

## Workflow
1. `docs/05_interpretation.md` 가이드 확인
2. `stages/04_analysis.md`에서 핵심 결과 확인
3. 결과의 도메인적 의미 해석
4. 한계점 및 후속 연구 방향 정리
5. 산출물 → `outputs/interpretation/run_XX/`에 저장
6. 아래 지식 업데이트
7. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 최종 리포트 작성 완료
- 한계점 및 후속 방향 명시됨
- 식별된 실패 원인이 다음 연구에 전달 가능한 형태로 정리됨

---

## 검증된 핵심 지식

### BioEval 해석 결과 (run_13~20, 2026-04-30 ~ 2026-05-04)

**H1: 순위 반전은 실재한다** (STRONG)
- MSE/R²와 BioEval-Dir은 모델을 독립적으로 순위 매긴다 (τ < 0.7)
- RPE1: τ(MSE, Dir_deg) = 0.389, 95% CI가 0 포함 → 통계적 독립
- Norman: mean_predictor가 MSE #1이나 Dir_deg #11 — Mean-Effect Trap 직접 증거
- DEG 비율이 낮을수록 Trap 심화 (Norman 1.53% < K562 2.38% < RPE1 6.50%)
- GEARS 포함 시 Norman만 τ=0.667 REVERSAL — 모델 다양성에 의존

**H2: BioEval은 domain-specific 이점** (MODERATE)
- Intra-DEG (DEG_auprc ↔ f1@50): 9/9 (100%), mean gap = +0.319
- Intra-magnitude (mag_rank ↔ gene-set): 6/6 (100%), mean gap = +0.303
- Cross-domain (Dir ↔ gene-set): 8/24 (33.3%) — 방향과 gene-set은 본질적 다른 신호
- MSE 자체가 domain-general predictor (ρ(-MSE, dir_discovery) = 0.88-0.96)
- BioEval 이점은 해석 가능성과 도메인 내 정밀도

**H3: 학습 모델 > baseline** (STRONG, qualified)
- Well-trained Ridge > baselines: 3 데이터셋 × 6 지표 ALL WIN
- GEARS(DL) < Ridge: 0/12 승 — model quality > model complexity
- GEARS vs mean_predictor: 11/12 승 — 순진한 baseline에는 승리
- Norman GEARS R²=-0.699, 128/283 유효 예측 — 'ctrl+gene' 형식 비호환
- BioEval이 품질 격차를 정확히 식별 — 판별력 검증

**핵심 인사이트**:
1. Mean-Effect Trap: MSE가 평균 근처 예측에 보상 → 방향 정보 상실
2. MSE의 이중성: domain-general predictor이나 방향 보장 불가
3. DL ≠ better: GEARS(GNN+attention)가 Ridge(선형)에 전패 — 복잡도가 품질을 보장하지 않음
4. 교세포 일관성: K562↔RPE1 모든 지표 τ > 0.78 — 현상이 체계적
5. DEG_auprc가 가장 견고한 H2 지표 — Bootstrap CI에서 모든 데이터셋 유의

**한계**:
- 4~9개 모델로 순위 해상도 제한 — GEARS 포함 시 K562/RPE1 τ=1.0
- GEARS 5 epochs + 최소 hyperparameter — 최대 성능 미대표 가능
- 단일 DL 모델 — CPA 등 추가 검증 필요
- 3개 인간 세포주 데이터셋 — 다른 생물종/조직 일반화 불확실

### 이전 연구에서 전달된 실패 원인
| 원인 | 단계 | 유형 | 전달 내용 |
|------|------|------|----------|
| Loss 불균형 (50:1) | Planning | prior | multi-objective loss는 반드시 분리 최적화 |
| Encoder 동결 | Planning | prior | pretrained component fine-tuning 필수 |
| 평가 오류 (leave-MoA-out 분류) | Framing | prior | 평가지표는 태스크 정의와 일치해야 함 |
| Surrogate R2 < 0 (135샘플) | Analysis | prior | GNN+XGBoost는 샘플 수 >= 500 필요 |
| Euler dFBA 불안정 | Analysis | prior | stiff ODE에는 implicit solver(BDF/Radau) 필수 |
| AL uncertainty 무효 | Analysis | prior | 모델 R2 > 0.3 이후에만 uncertainty 기반 AL 사용 |

---

## 클레임 카드 (claim_card.yaml)

| ID | 가설 | 클레임 | 강도 |
|----|------|--------|------|
| C1 | H1 | MSE/R²와 BioEval-Dir은 섭동 예측 모델을 통계적으로 독립적으로 순위 매긴다. 순위 불일치는 모델 다양성이 확보될 때 가장 뚜렷하며, DEG 비율이 낮은 데이터셋에서 더 심화한다. | STRONG |
| C2 | H2 | BioEval 지표는 intra-domain에서 MSE보다 downstream 생물학적 유용성을 더 잘 예측한다. Cross-domain에서는 MSE가 더 나은 predictor이다. BioEval의 가치는 domain-specific 정밀도와 해석 가능성에 있다. | MODERATE |
| C3 | H3 | 잘 학습된 모델(Ridge)은 모든 BioEval 지표에서 baseline을 압도한다. DL 복잡도 자체가 우위를 보장하지 않는다 — GEARS(DL)는 Ridge(선형)에 전패. 이 결과 자체가 BioEval의 판별력을 검증한다. | STRONG (qualified) |

---

## Run 이력 (세부 내용은 outputs/interpretation/run_XX/ 참조)
- **run_01** (2026-05-03): BioEval 전체 해석 리포트 — H1 STRONG, H2 MODERATE (domain-specific), H3 STRONG (qualified). GEARS < Ridge 확인. Mean-Effect Trap 기작 해석. MSE 이중성 분석.
- **claim_card** (2026-05-04): C1(H1 STRONG), C2(H2 MODERATE domain-specific), C3(H3 STRONG qualified) 작성 완료
