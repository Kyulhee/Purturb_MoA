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

### 이전 방향(NAP)에서 검증된 지식 (참고용)
- **XGBoost-only R2=0.91** (textbook, 137차원 knockout mask) — FBA는 근본적으로 tabular problem
- **NAP 가치 조건**: textbook 0/6, multi-species 2/6, FlowGAT-style 2/6 — novelty 부족으로 방향 전환
- **AL 실패**: input-space AL R2=0.69 vs random R2=0.81 — ensemble uncertainty가 FBA surrogate에 부적합
- **dFBA 안정성**: BDF 정확, Euler 사용 금지 — 후속 연구 참고
- **GPU 벤치마크**: textbook CPU 우세, iJO1366 GPU 1.3x — 소규모 그래프는 CPU 효율

### Perturb-seq 미해결 문제 (run_04에서 검증)
- **Problem 1**: 평가 지표 붕괴 — Islander(MLP)가 SOTA 능가하나 생물학 왜곡 (Wang et al., 2026)
- **Problem 2**: 조합 폭발 — 단일 섭동만으로 조합 예측하는 방법 없음
- **Problem 3**: 교차 세포 유형 전이 — GEARS 명시적 미지원, Cell-JEPA 효과 크기 불가
- **Problem 4**: 근본적 난이도 이론 부재 — 정보이론적 가변성 한계 없음
- **Problem 5**: 측정-필요 불일치 — immortalized KO vs primary drug response

### 교차 도메인 탐색 결과 (run_04에서 검증)
- OT 도메인: 포화 — avoid
- MapPFN(arXiv:2601.21092): 인컨텍스트 학습으로 제로샷 가능, 2026년 1월 — 최신
- FCR+ICM 결합: 미탐색, 참신성 VERY HIGH
- 조합 일반화: 원칙적 프레임워크 부재, 참신성 VERY HIGH

---

## 다음 단계

1. ~~Step 0-3: 기존 분석(NAP) 재현~~ 완료 (novelty 부족으로 방향 전환)
2. **Step 1**: 데이터 확보 및 전처리 (Norman 2019, Replogle 2022)
3. **Step 2**: FCR 인코더 구현 + z_x/z_t/z_tx 분해 검증
4. **Step 3**: ICM 정규화 구현 + z_tx 불변성 검증 (RQ1)
5. **Step 4**: 경로 모듈 분해 + 조합 함수 구현 (RQ2)
6. **Step 5**: 교차 세포 유형 제로샷 전이 실험 (RQ3)
7. **Step 6**: 소거 실험 매트릭스 (6개 구성)
8. **Step 7**: 결과 분석 + 논문 초안

---

## Run 이력 (세부 내용은 outputs/analysis/run_XX/ 참조)
- **run_01** (2026-04-26): Module A/B/C 구현 + E2E 파이프라인. XGBoost R2=0.91, NAP 가치 조건 매핑, dFBA 안정성, GPU 벤치마크
- **run_02** (2026-04-27): 문헌 심층 리뷰 11편. FlowGAT, scFEA, SARTRE, Xu 2019, homophily
- **run_03** (2026-04-27): Input-space AL 실험. AL R2=0.69 vs Random R2=0.81 — AL 실패
- **run_04** (2026-04-27): NAP transfer v2 실험(불완전) + Perturb-seq 미해결 문제 43편 리뷰 + 교차 도메인 10개 스캔 → 방향 전환 계기
