# Stage 03 — Planning

## Loopback 기록
- **2026-04-30**: Framing → Planning. BioEval 프로젝트(섭동 예측 평가 지표) 실험 설계. 4-Phase: 데이터 확보→지표 구현→순위 반전 분석→downstream 상관 분석. 상세: `outputs/planning/run_06/`

## Workflow
1. `docs/03_planning.md` 가이드 확인
2. `stages/02_framing.md`에서 베이스라인/타겟 확인
3. 실험 설계 작성 (모델, 피처, 하이퍼파라미터, 평가 전략)
4. **사용자 컨펌 획득** (타겟 성능 + 실험 설계)
5. 산출물 → `outputs/planning/run_XX/`에 저장
6. 아래 지식 업데이트 (과거 run에서 검증된 인사이트만 통합)
7. `objects/current/experiment_contract.yaml` 업데이트
8. `objects/current/evaluation_validity_card.yaml` 업데이트
9. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 실험 설계서 작성 완료
- 사용자 컨펌 획득 ✅ confirmed by user 2026-04-30

---

## 연구 질문 (from stages/02)

**"섭동 예측 평가 지표가 생물학적 충실도를 측정하는가? 생물학적 충실도를 측정하는 지표를 설계하고, 이 지표 하에서 Ahlmann-Eltze의 'DL ≤ baseline' 위기가 해소되는가?"**

### 하위 질문
1. **RQ1 (지표-생물학 상관 진단)**: 기존 지표 중 어떤 것이 생물학적 유용성과 상관하는가?
2. **RQ2 (지표 설계)**: BioEval 지표가 기존 지표보다 downstream 생물학적 유용성을 더 잘 예측하는가?
3. **RQ3 (베이스라인 위기 해소)**: BioEval 하에서 DL > baseline 체계가 존재하는가?

---

## 1. 방법론 아키텍처 (4-Phase)

```
Phase 1: 데이터 확보 + 모델 예측 수집
Phase 2: BioEval 지표 구현 + 기존 지표 계산
Phase 3: 지표-순위 반전 분석 (RQ1 + RQ3)
Phase 4: 지표-downstream 과업 상관 분석 (RQ2)
```

## 2. 설계 결정과 근거

| 설계 결정 | 근거 | 대안 |
|-----------|------|------|
| Kendall τ를 RQ3 주 지표로 | 비모수 순위 통계, 순환 위험 없음 | AUROC(순환 위험, run_09 교훈) |
| BioEval-Dir: 유전자×섭동 수준 분해 | 기존 지표는 섭동 수준만(Shesha) 또는 유전자 수준만(AUPRC) | 전체 평균(분해 정보 손실) |
| 다수 모델 비교 (5개) | 순위 분석에 최소 4-5개 모델 필요 | 단일 모델(순위 분석 불가) |
| Ahlmann-Eltze 예측 확보 시도 → 실패 시 직접 학습 | 재현성 우선, 직접 학습은 대안 | 직접 학습만(원 논문과 비교 불가) |
| 둘 다 발표 가능 설계 | τ < 0.5(반전) 또는 τ > 0.7(유지) 모두 중요 결과 | 단방향 성공만 목표(위험) |
| S1-S7 민감도/소거 실험 | BioEval 구현 민감성 방어 | 임계값 고정(재현성 위험) |

## 3. 평가 타겟

| RQ | 지표 | Baseline | 타겟 | 데이터 |
|----|------|----------|------|--------|
| RQ1 | 지표-AL 상관 (Spearman) | 0 (MSE-AL 상관) | > 0.5 | Replogle |
| RQ1 | 지표-DEG 상관 (Spearman) | AUPRC 기준치 | BioEval > MSE by ≥0.1 | Replogle |
| RQ2 | BioEval-Dir 방향 정확도 | 0.5 (우연) | > 0.7 | Replogle, Norman |
| RQ2 | BioEval-Cal 보정 기울기 | — | 0.8-1.2 | Replogle |
| RQ3 | Kendall τ (MSE vs BioEval) | 1.0 (동일) | < 0.5 또는 > 0.7 | Replogle, Norman |
| RQ3 | Synergistic GI DL>baseline | 0% (MSE 하) | > 30% (BioEval 하) | Norman |

## 4. 소거/민감도 실험

| ID | 변인 | 목적 |
|----|------|------|
| S1 | DEG 임계값 스윕 (0.1, 0.25, 0.5, 1.0) | BioEval-Dir 견고성 |
| S2 | BioEval-Composite 가중치 변동 | 통합 지표 견고성 |
| S3 | 방향 평가: 전체 vs DEG-only | 유전자 선택 효과 |
| S4 | 보정 분석 분해 수준 | 분해 수준 효과 |
| S5 | 데이터셋 교차 Replogle↔Norman | 일반성 |
| S6 | Baseline 모델 추가/제거 | 순위 안정성 |
| S7 | AUPRC vs BioEval-DEG 직접 비교 | 기존 DEG 지표와 차이 |

## 5. 데이터

| 데이터셋 | 용도 | 비고 |
|----------|------|------|
| Replogle 2022 | RQ1-3 주 평가 | K562+RPE1, 848 공유 섭동 |
| Norman 2019 | RQ1-3 조합 평가 | 128 double-KO, GI ground-truth |
| PBMC (Zhu 2025) | RQ1 AUPRC 비교 | 7 cell types, IFN-γ |
| PORTAL 2026 | RQ3 대규모 검증 (선택) | 665K pairwise |

## 6. 핵심 리스크

| 리스크 | 확률 | 영향 | 대응 |
|--------|------|------|------|
| Ahlmann-Eltze 예측 확보 불가 | 중간 | 높음 | 직접 모델 학습으로 대체 |
| 순위 반전이 안 일어남 | 중간 | 낮음 | 부정 결과도 발표 가능 |
| BioEval 구현 민감성 | 중간 | 중간 | S1-S4 민감도 실험 |
| DEG 임계값이 결과 지배 | 중간 | 중간 | S1 임계값 스윕 |

## 7. 이전 실패에서의 설계 원칙

| 교훈 | 출처 | 반영 |
|------|------|------|
| 지표 선택이 결론 변경 | run_12 | RQ의 직접적 동기 |
| 경쟁자 조기 확인 | FCR-ICM | Framing에서 5개 부분 경쟁자 확인 |
| 잠재공간-유전자공간 갭 | run_05 | 모든 지표 유전자 공간에서 계산 |
| 소거실험 필수 | run_09 | S1-S7 민감도 매트릭스 |
| AUROC=1.0 동어반복 | run_09 | Kendall τ(순위 상관)로 순환 방지 |

---

## Run 이력 (세부 내용은 outputs/planning/run_XX/ 참조)
- **run_06** (2026-04-30): BioEval 실험 설계. 4-Phase, 3 RQ, 7 소거/민감도 실험. Ahlmann-Eltze 예측 확보 시도 + 직접 학습 대안
- run_01-05 (2026-04-25~30): [FCR-ICM] 이전 프로젝트 실험 설계
