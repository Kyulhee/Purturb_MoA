# Stage 04 — Analysis

## Workflow
1. `docs/04_analysis.md` 가이드 확인
2. `stages/03_planning.md`에서 실험 설계 및 기준 확인
3. 데이터 전처리, 모델 학습, 평가 수행
4. 결과를 Planning 타겟과 비교
5. 산출물 → `outputs/analysis/run_XX/`에 저장
6. 아래 지식 업데이트
7. CLAUDE.md의 `current_stage` 업데이트

## Done when
- Planning 타겟 성능 달성 또는 미달 시 사용자 보고 후 방향 결정

---

## 검증된 핵심 지식

### H1 인과 불변성 가설 — 실험 결과 (run_04-05)

**RQ1: ICM이 z_tx를 cell type 불변으로 만드는가** — 합성 통과
- FCR (ICM 없음): 교세포 상관계수 0.505
- FCR + ICM: 교세포 상관계수 0.971
- ICM MMD 정규화가 10개 교란 모두에서 개선 (범위: +0.21 ~ +0.74)

**RQ2: 단일-KO z_tx로 조합적 예측 가능한가** — 합성(잠재공간) 실패, 합성(유전자공간) 통과, 실제 통과
- 합성 잠재공간: best_corr=0.29, best_R2=-1.63 — 인코더 비선형 변환으로 조합 구조 파괴
- 합성 유전자공간: **R2=0.88, corr=0.94** — 디코더가 인코더 비선형성 보상 (run_05로 검증)
- Norman 2019 실제: FCR best_corr=0.955, best_R2=0.881 / FCR+ICM best_corr=0.951, best_R2=0.870
- **소거실험(run_05)**: 조합 일관성 손실이 RQ2-cross 0.20→0.79 개선. 전체 모델(config 6) RQ1=0.99, RQ3=0.99
- **실제데이터 소거실험(run_06)**: Norman 128 double-KO 쌍에서 모든 config R2=0.86-0.89로 유사. comp loss가 오히려 성능 저하(0.88→0.86). baseline FCR이 이미 조합성에 충분

**RQ3: ICM이 zero-shot 교세포 전이를 가능하게 하는가** — 합성 통과
- FCR (ICM 없음): 전이 상관계수 0.508, 코사인 0.476
- FCR + ICM: 전이 상관계수 0.960, 코사인 0.956

**Replogle 2022**: K562 단일 세포유형만 확보 → RQ1/RQ3 다세포유형 검증 불가

### 핵심 인사이트
1. **ICM은 불변성에 확실히 유효** (RQ1, RQ3). MMD 정규화가 z_tx 교세포 정렬에 강력
2. **RQ2 합성-실제 갭 해명됨**: 조합성은 유전자 공간에서 평가해야 함. 잠재 공간 R2=0.05 vs 유전자 공간 R2=0.88. 디코더가 인코더 비선형성 보상
3. **조합 일관성 손실이 핵심 구성요소**: ICM만으로 RQ2 개선 안 됨. comp_loss가 RQ2-cross를 0.20→0.79로 향상
4. **ICM이 인코더를 더 선형적으로 만듦**: linear R2 0.69→0.87. ICM의 분포 정렬이 부차적으로 조합 구조 보존에 기여
5. **실제데이터에서 comp loss 불필요**: Norman 소거실험에서 baseline FCR이 이미 R2=0.88. comp loss는 과도한 제약으로 오히려 저하
6. **가법성이 우위**: 128 double-KO 중 110쌍(86%)에서 가법 조합이 승법보다 우수

### 이전 방향(NAP)에서 검증된 지식 (참고용)
- **XGBoost-only R2=0.91** — FBA는 근본적으로 tabular problem
- **GNN 임베딩 중복**: 정적 그래프에서는 knockout mask가 충분 통계량
- **AL 실패**: FBA가 싸고 입력 차원이 낮아 AL 이점 없음

---

## 현재 진행 상태

| 단계 | 상태 | 비고 |
|------|------|------|
| Step 1: 데이터 확보 | 완료 | Norman 2019 (89K cells), Replogle K562 (162K cells) |
| Step 2: FCR 인코더 구현 | Phase 1 완료 | VAE, z_dim=8, z_x/z_t/z_tx 분해 |
| Step 3: ICM 정규화 (RQ1) | 합성 통과 | 교세포 corr 0.505→0.971 |
| Step 4: 조합성 (RQ2) | 합성 실패, 실제 통과 | Norman best_corr=0.955, best_R2=0.881 |
| Step 5: 교세포 전이 (RQ3) | 합성 통과 | 교세포 corr 0.508→0.960 |
| Step 6: 소거 실험 | 완료 (run_05) | 6구성 소거 + 인코더 비선형성 + 갭 해명 |
| Step 7: 논문 초안 | 미실행 | |

## 미해결 과제

1. **RQ1/RQ3 실제데이터 검증**: 다세포유형 Perturb-seq 데이터셋 필요 (Replogle RPE1 로딩 실패)
2. ~~**RQ2 합성-실제 갭 해명**~~: run_05에서 해명 — 유전자 공간 평가로 전환 시 디코더 보상으로 R2=0.88 달성
3. ~~**소거 실험**~~: run_05 완료 — ICM 핵심, comp_loss 필수, linear head 선택적
4. **실제데이터 소거실험**: 완료 (run_06) — baseline FCR이 이미 R2=0.88, comp loss 불필요 확인
5. **논문 초안 작성**: Stage 05로 이관

---

## Run 이력 (세부 내용은 outputs/analysis/run_XX/ 참조)
- **run_01** (2026-04-26): NAP E2E 파이프라인. XGBoost R2=0.91, GNN 중복성 확인
- **run_02** (2026-04-27): GNN vs tabular 문헌 심층 리뷰 11편
- **run_03** (2026-04-27): Input-space AL 실험. AL R2=0.56 vs Random R2=0.68 — AL 실패
- **run_04** (2026-04-27): Perturb-seq 미해결 문제 리뷰 + 교차 도메인 10개 스캔 → 방향 전환 + Phase 1 합성검증(RQ1/RQ3 통과, RQ2 실패) + Phase 2 Norman 실제데이터(RQ2 통과: best_corr=0.955, best_R2=0.881) + Replogle RPE1 로딩 실패
- **run_05** (2026-04-27): 소거실험(6구성) + 인코더 비선형성 측정 + RQ2 갭 해명(잠재공간 R2=0.05 vs 유전자공간 R2=0.88)
- **run_06** (2026-04-27): Norman 실제데이터 소거실험 — 모든 config R2=0.86-0.89, comp loss 오히려 저하, 가법성 86% 우위
