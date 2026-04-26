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

### 치명적 실패 원인 (run_01, 재발 방지)
**3중 실패: Framing + Planning + Data** — 이 중 하나라도 해결 안 되면 재시도 불가

1. **Framing**: leave-MoA-out 분류 accuracy → 구조적 0% 보장 → 클러스터링 품질로만 평가
2. **Planning**: L_recon 99.5% gradient → MoA 분리 무시 → effective gradient 기반 weight 필수
3. **Data**: 합성 데이터만 사용 → 신호 약함 → **실제 sci-Plex 데이터 필수**

### 학습 역설 (run_01에서 관찰)
- accuracy가 학습 진행에 따라 하락 (84%→76%)
- 재구성 개선 ↔ MoA 분리 악화 → loss 불균형의 직접적 증거

### 재시도 전제조건 (모두 충족 필요)
1. 실제 sci-Plex 데이터 사용
2. Loss weight effective gradient 기반 설계
3. Drug encoder 동결 해제 또는 LR 차등 적용
4. Leave-MoA-out = 클러스터링 품질 평가

### V2 결과 (참고용, 합성 데이터)
| 분할 | Cls Top-1 | Metric Top-1 | Zero-Shot |
|------|-----------|-------------|-----------|
| random | 0.76 | 0.55 | 1.00 |
| leave_compound_out | 0.13 | 0.09 | 0.09 |
| leave_moa_out | 0.00 | 0.00 | 0.00 |

---

## Run 이력 (세부 내용은 outputs/analysis/run_XX/ 참조)
- run_01: 합성 데이터 Perturb-seq MOA. LCO 9%, LMO 0%. 3중 실패(Framing+Planning+Data) → Literature Review로 루프백
