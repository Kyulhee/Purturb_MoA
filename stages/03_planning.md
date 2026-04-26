# Stage 03 — Planning

## Workflow
1. `docs/03_planning.md` 가이드 확인
2. `stages/02_framing.md`에서 베이스라인/타겟 확인
3. 실험 설계 작성 (모델, 피처, 하이퍼파라미터, 평가 전략)
4. **사용자 컨펌 획득** (타겟 성능 + 실험 설계)
5. 산출물 → `outputs/planning/run_XX/`에 저장
6. 아래 지식 업데이트 (과거 run에서 검증된 인사이트만 통합)
7. CLAUDE.md의 `current_stage` 업데이트

## Done when
- 실험 설계서 작성 완료
- 사용자 컨펌 획득 (`confirmed by user [날짜]`)

---

## 검증된 핵심 지식 (모든 run에서 통합)

### 1. 치명적 설계 결함 3가지 (run_01 실패 원인, 반복 금지)
| 결함 | 교훈 | 올바른 설계 |
|------|------|------------|
| Loss 불균형 (50:1) | effective gradient 크기로 weight 설계 | L_recon weight 0.01~0.05 또는 MoA loss 10~50배 증폭 |
| Drug encoder 동결 | 사전학습 후에도 fine-tuning 필요 | LR 차등 적용(1:10~1:100) 또는 동결 해제 |
| leave-MoA-out 분류 평가 | 학습되지 않은 라벨 분류 불가 | 클러스터링 품질(silhouette, ARI, NMI)로 평가 |

### 2. 아키텍처 유효 성분
- rFCFP 인코딩, metric prototype head, 대조 손실 아이디어: 유효
- z_perturb.detach(): 재구성/MoA gradient 분리에 필요
- 구현 디테일(loss 비율, 동결, 평가)에서 실패했을 뿐 구조는 타당

### 3. Phase 3-4 통합 아키텍처

```
Module A: FBA Ground Truth (COBRApy 병렬, 10K조합 ~10-30분)
    → Module B: GNN+XGBoost Surrogate (HGTConv 3층 + XGBoost)
    → Module C: Active Learning Loop (UCB/EI, FBA 호출 70-90% 감소)
    → Module D: COMETS dFBA (Java 필수, dfba-python 대안)
    → Module E: NSGA-II (pymoo, 접종 비율 최적화, 3목적)
    → Module F: TOPSIS/Pareto (Entropy weight + 민감도 분석)
```

**의존 구조**: A→B→C (Phase 3), D→E→F (Phase 4), C↔D (Active Learning 양방향)

### 4. Phase 실현성 평가
- **Phase 3**: MEDIUM-HIGH (GNN 임베딩 동적 특성 포착, 일반화 리스크)
- **Phase 4**: MEDIUM (FLYCOP 불가→COMETS 대체, dFBA 수치 안정성, NSGA-II 계산 비용)

### 5. FLYCOP 대체 전략 (저장소 삭제 확인)
- FLYCOP: 저장소 삭제, 사용 불가. 원논문 전문도 접근 불가
- 대체: COMETS v2.12.4 (dFBA+공간), dfba-python (비공간 순수 Python)
- FLYCOP의 fuzzy 다목적 평가 → TOPSIS + Entropy weight로 객관화 대체

### 6. COMETS 설치 검증 결과
- cometspy v0.6.3: pip 설치 성공, Python 인터페이스 정상
- **COMETS Java 코어 필수**: COMETS_HOME 환경변수 설정 필요
- 대안: Docker 컨테이너 활용 또는 dfba-python (순수 Python, Java 불필요)

### 7. 오픈소스 도구
| 도구 | 버전 | 용도 | URL |
|------|------|------|-----|
| COBRApy | 0.31.1 | FBA | https://github.com/opencobra/cobrapy |
| PyTorch Geometric | latest | 이종 GNN | https://github.com/pyg-team/pytorch_geometric |
| XGBoost | latest | 회귀 | https://github.com/dmlc/xgboost |
| COMETS | 2.12.4 | dFBA | https://github.com/segrelab/COMETS |
| cometspy | 0.6.3 | COMETS Python | https://github.com/segrelab/cometspy |
| pymoo | latest | NSGA-II | https://github.com/anyoptimization/pymoo |

### 8. 리소스 관리 규칙
- RAM < 4GB, Disk < 10GB 시 작업 중단 + 사용자 보고
- Multi-Agent ≤ 3개 동시 (RAM 40GB 기준)
- GNN/FBA 메모리 집약 에이전트는 1개만

---

## 다음 단계 (미해결)
1. COMETS Java 코어 설치 또는 dfba-python으로 dFBA 실행 검증
2. GNN 대사 네트워크 변환 프로토타입 완료 (백그라운드 진행 중)
3. COBRApy 병렬 FBA 벤치마크 완료 (백그라운드 진행 중)
4. 실제 sci-Plex 데이터 사용 (합성 데이터 한계 확인됨)
5. Loss weight 재설계: effective gradient 기반
6. 사용자 컨펌 절차 의무화

---

## Run 이력 (세부 내용은 outputs/planning/run_XX/ 참조)
- run_01: 임계값 분류+MoA 계획 → 3가지 치명적 결함으로 Analysis 전면 실패
- run_02: Phase 3-4 기술실현성 분석 → FLYCOP 삭제 확인, COMETS 대체
- run_03: 통합 아키텍처 설계 + 6모듈 정의 + 리소스 규칙 + FLYCOP/COMETS 분석
