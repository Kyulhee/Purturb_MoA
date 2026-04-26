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

### 이전 연구에서 전달된 실패 원인
| 원인 | 단계 | 전달 내용 |
|------|------|----------|
| Loss 불균형 (50:1) | Planning | multi-objective loss는 반드시 분리 최적화 |
| Encoder 동결 | Planning | pretrained component fine-tuning 필수 |
| 평가 오류 (leave-MoA-out 분류) | Framing | 평가지표는 태스크 정의와 일치해야 함 |
| Surrogate R2 < 0 (135샘플) | Analysis | GNN+XGBoost는 샘플 수 >= 500 필요 |
| Euler dFBA 불안정 | Analysis | stiff ODE에는 implicit solver(BDF/Radau) 필수 |
| AL uncertainty 무효 | Analysis | 모델 R2 > 0.3 이후에만 uncertainty 기반 AL 사용 |

---

## Run 이력 (세부 내용은 outputs/interpretation/run_XX/ 참조)
