# Stage 02 — Framing

## Loopback 기록
- **2026-04-30**: 신규 프로젝트 전환. FCR-ICM 프로젝트(run_12 PARTIAL)에서 BioEval(섭동 예측 평가 지표)로 완전 전환. 상세: `outputs/framing/run_06/`

## Workflow
1. `docs/02_framing.md` 가이드 확인
2. `stages/01_literature_review.md`에서 인사이트 확인
3. 연구 질문, 베이스라인, 타겟 성능 수치 정의
4. 산출물 → `outputs/framing/run_XX/`에 저장
5. 아래 지식 업데이트 (검증된 인사이트만 통합)
6. `objects/current/novelty_ledger.yaml` 업데이트 (신규성 4계층, 기여 분해, 경쟁 밀도, 허용 클레임 강도)
7. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 연구 질문이 단일 문장으로 기술 가능
- 베이스라인 수치가 문헌 근거와 함께 명시됨
- 타겟 성능 수치가 베이스라인 기반으로 설정됨
- 평가 지표 및 데이터셋 확정

---

## 연구 질문

**"섭동 예측 평가 지표가 생물학적 충실도를 측정하는가? 생물학적 충실도를 측정하는 지표를 설계하고, 이 지표 하에서 Ahlmann-Eltze의 'DL ≤ baseline' 위기가 해소되는가?"**

### 하위 질문

1. **RQ1 (지표-생물학 상관 진단)**: 기존 평가 지표(MSE, R², Pearson, DEG overlap, PDS) 중 어떤 것이 생물학적 유용성(DEG 회복, 방향 정확도, downstream 과업 성과)과 상관하는가?
2. **RQ2 (지표 설계)**: 생물학적으로 기반하고, 유전자 수준 분해능을 가지며, 방향을 인식하고, 구현에 견고한 평가 지표를 설계할 수 있는가? 이 지표가 기존 지표보다 downstream 생물학적 유용성을 더 잘 예측하는가?
3. **RQ3 (베이스라인 위기 해소)**: 생물학적 충실도 지표 하에서 DL 모델이 simple baseline을 능가하는 체계(synergistic GI, 교세포 전이, 고다기능 섭동)가 존재하는가? 아니면 베이스라인 위기가 지표와 무관하게 실재하는가?

### 도출 근거

**해결하는 Gap (문헌 리뷰 run_07 Gap 2)**:
- SCALE (2026): MSE가 mean-effect trap 유발 — 지표가 모델 행동을 왜곡
- Ahlmann-Eltze (2025): DL ≤ linear baseline across 7+ benchmarks — 위기 원인 불명
- Shesha (2026): magnitude와 stability가 0.75-0.97 상관하나 불일치 사례에서 생물학 노출
- 우리 run_12: prod rho=0.437 PASS vs A7 rho=0.326 PARTIAL — 지표가 결론 변경

**핵심 차별화**: Wei et al. (2026)은 "어떤 모델이 좋은가?"를 물음. 우리는 **"지표를 바꾸면 모델 순위가 어떻게 변하는가?"**를 물음. 이 메타-평가 질문은 누구도 묻지 않음.

### 평가 전략

| RQ | 지표 | 베이스라인 | 타겟 | 데이터 |
|----|------|----------|------|--------|
| RQ1 | 지표-AL 상관 (Spearman) | 0 (MSE-AL 상관) | > 0.5 | Replogle |
| RQ1 | 지표-DEG 상관 (Spearman) | AUPRC 기준치 | BioEval > MSE by ≥0.1 | Replogle, PBMC |
| RQ2 | BioEval-Dir 방향 정확도 | 0.5 (우연) | > 0.7 | Replogle, Norman |
| RQ2 | BioEval-Cal 보정 기울기 | — | 0.8-1.2 범위 | Replogle |
| RQ3 | Kendall τ (MSE vs BioEval 순위) | 1.0 (동일) | < 0.5 (반전) 또는 > 0.7 (유지) | Replogle, Norman |
| RQ3 | Synergistic GI에서 DL>baseline 비율 | 0% (MSE 하) | > 30% (BioEval 하) | Norman |

### 베이스라인 계층

1. **우연 수준**: 랜덤 예측 (MSE 최대, BioEval-Dir=0.5, AUPRC=π₀)
2. **Mean predictor**: 관측치 평균 (Ahlmann-Eltze 베이스라인)
3. **Additive linear**: Y = GW^T P + b (Ahlmann-Eltze 최우수)
4. **scGPT + linear**: 선형 모델 + 사전학습 임베딩 (Ahlmann-Eltze에서 최고)
5. **CPA**: 조합 오토인코더 (표준 비교)
6. **GEARS**: GNN + GRN (조합 예측 표준)
7. **SCALE**: LLaMA 기반 파운데이션 모델 (최신 FM)

### 데이터

| 데이터셋 | 용도 | 비고 |
|----------|------|------|
| Replogle 2022 | RQ1-3 주 평가 | K562+RPE1, 848 공유 섭동, 교세포 |
| Norman 2019 | RQ1-3 조합 평가 | 128 double-KO, GI ground-truth |
| Ahlmann-Eltze 벤치마크 | RQ3 직접 재현 | 7+ 벤치마크, 모델 예측 확보 필요 |
| PBMC (Zhu 2025) | RQ1 AUPRC 비교 | 7 cell types, IFN-γ 자극 |
| PORTAL 2026 | RQ3 대규모 검증 | 665K pairwise (선택) |

### 신규성 검사 결과

| 검사 | 결과 | 상세 |
|------|------|------|
| 경쟁 밀도 | 5개 부분 경쟁자 | Wei(벤치마크), Zhu(AUPRC), SCALE(Cell-Eval), Shesha(안정성), Csendes(FM≤mean). 직접 동일 접근(지표 설계+순위 반전+downstream 상관)은 0개 |
| 기여 분해 | 3/6 성분 독립 신규 | BioEval-Dir(유전자×섭동 분해), BioEval-Cal(보정 분석), 지표-순위 반전 정량화, 지표-downstream 과업 상관 |
| 실패 모드 | 2개 ⚠ | 포화 시장(5개 부분 경쟁자), 평가 함정(BioEval 자체 구현 민감성) |
| 둘 다 발표 가능 | 예 | 순위 반전→"지표 아티팩트"(높은 임팩트), 순위 유지→"위기 실재"(중간 임팩트) |

### 핵심 위험

1. **Wei et al. (2026)과 중복**: Wei는 기존 지표로 벤치마크, 우리는 새 지표 설계+순위 반전 분석. 본질적으로 다른 과업
2. **BioEval 구현 민감성**: SCALE이 Cell-Eval에서 지적한 문제. 임계값 스윕으로 견고성 테스트 필요
3. **순위 반전이 안 일어남**: 부정 결과도 발표 가능. "위기 실재" = 중요 발견
4. **Ahlmann-Eltze 예측 확보 불가**: 코드 공개 확인 필요. 미공개 시 Replogle+Norman으로 재현

---

## Run 이력 (세부 내용은 outputs/framing/run_XX/ 참조)
- **run_06** (2026-04-30): BioEval 신규 프로젝트 프레이밍. 3 RQ(지표-생물학 상관, 지표 설계, 베이스라인 위기 해소). 핵심 차별화: 지표-순위 반전 분석(Kendall τ). 5개 부분 경쟁자, 3/6 독립 신규 성분. FCR-ICM에서 전환
- **run_05** (2026-04-30): [FCR-ICM] Analysis 루프백. RQ3를 "ICM 기반 전이"→"잔차 순위 기반 전이"로 재구성
- **run_04** (2026-04-29): [FCR-ICM] Analysis 루프백. RQ3 동어반복 해소
- **run_03** (2026-04-28): [FCR-ICM] A+B 결합 프레이밍
- **run_01-02** (2026-04-25~27): [FCR-ICM] 초기 프레이밍
