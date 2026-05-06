# Analysis — Guide

## Pre-Flight Resource Check (필수)
분석/데이터 처리 작업 시작 전 반드시 시스템 리소스를 확인할 것:

```
python -c "import psutil; m=psutil.virtual_memory(); d=psutil.disk_usage('C:/'); print(f'RAM: {m.available/(1024**3):.1f}/{m.total/(1024**3):.1f}GB ({m.percent}%) | Disk: {d.free/(1024**3):.1f}/{d.total/(1024**3):.1f}GB ({d.percent}%) | CPU: {psutil.cpu_count(logical=False)}P/{psutil.cpu_count(logical=True)}L')"
```

| 리소스 | 정상 | 경고 | 중단 |
|--------|------|------|------|
| Available RAM | >8GB | 4-8GB (가벼운 작업만) | <4GB (사용자 보고) |
| Available Disk | >20GB | 10-20GB (필요시 정리) | <10GB (사용자 보고) |
| CPU usage | <70% | 70-90% (경량 작업) | >90% (에이전트 축소) |

**Multi-Agent 실행 시 주의**:
- 동시 백그라운드 에이전트 ≤ 3개 (RAM 40GB 기준)
- GNN/FBA 등 메모리 집약적 에이전트는 1개만 실행
- OOM 방지: 대규모 데이터셋은 chunk 단위로 처리

## Implementation Process (구현 프로세스 지침)

모든 분석 작업은 아래 6단계 순서로 진행한다. 각 Step 완료 시 반드시
`outputs/analysis/run_XX/`에 체크포인트 파일(중간 결과, 상태 JSON)을 저장하여
LLM 세션 컨텍스트가 재설되더라도 이어서 작업 가능하도록 한다.

### Step 0: 데이터 접근성 검증 + 프로파일링
- 실제 데이터(GEO 다운로드, SBML 모델 파일 등) 접근 가능한지 선행 확인
- 다운로드 타임아웃, 포맷 오류, 라이선스 제약 등을 이 단계에서 발견
- 합성/샘플 데이터로 대체 가능한지, 불가능하면 Step 진행 중단
- **산출물**: 데이터 접근성 보고서 (`data_access_report.md`)

### Step 1: 기존 분석 재현
- 기존 베이스라인 코드를 그대로 실행하여 결과 재현
- 재현 실패 시 원인(환경 차이, 누락 의존성, 데이터 불일치) 기록
- **산출물**: 재현 결과 + 환경 로그 (`reproduction_results.md`)

### Step 2: 결과 재현 확인 (수치 기준 사전 정의)
- Step 1 시작 전 재현 성공 기준을 명시적으로 정의
  - 예: "베이스라인 논문의 Top-1 accuracy ±5% 이내"
  - 예: "FBA growth rate 오차 < 1%"
- 기준 미달 시 원인 분석 후 Step 1 재시도 또는 기준 조정 (사용자 승인)
- **산출물**: 재현 확인 보고서 (`reproduction_verification.md`)

### Step 3: 스몰 서브셋으로 모듈별 구현-검증
- 전체 파이프라인을 한 번에 구현하지 않고, 모듈 단위로 구현-검증 사이클 반복
- 모듈 순서는 Planning(stages/03)에서 정의한 파이프라인 아키텍처를 따른다
- 예시: A(데이터 생성) → B(모델 학습) → C(탐색 전략) → D(시뮬레이션) → E(최적화) → F(의사결정)
- 각 모듈에서 입력/출력 인터페이스(데이터 포맷, 스키마)를 먼저 확정
- 각 모듈 완료 시 단위 테스트 + 체크포인트 저장
- **산출물**: 모듈별 코드 + 단위 테스트 + 인터페이스 문서

### Step 4: 베이스라인 비교 + 모듈별 Ablation
- "베이스라인 vs 우리 모델" 전체 비교만으로는 어떤 변경이 효과적이었는지 파악 불가
- Ablation은 각 모듈의 핵심 설계 결정(design decision) 단위로 변형(variant)을 만들어 비교
- 변형은 Planning 단계에서 미리 정의하되, Analysis 중 새로운 비교 필요성이 발견되면 추가 가능
- Ablation 결과는 테이블 형태로 정리하여 해석력 확보
- **산출물**: ablation 결과 테이블 + 해석 보고서

### Step 5: 실제 규모 데이터로 벤치마크
- Step 0에서 검증된 실제 데이터로 전체 파이프라인 실행
- 스몰 서브셋에서 확인된 인터페이스/파라미터 그대로 사용
- 성능(정확도, 실행 시간, 메모리) 벤치마크 수집
- **산출물**: 최종 벤치마크 보고서 (`full_scale_benchmark.md`)

### 체크포인트 규칙
- 각 Step 완료 시 `outputs/analysis/run_XX/checkpoint_stepN.json` 저장
- 체크포인트 내용: Step 번호, 완료 시간, 핵심 결과 요약, 다음 Step 전제조건
- 다음 세션에서는 `stages/04_analysis.md` + 최신 체크포인트로 작업 재개

---

## Checklist
- [ ] **시스템 리소스 체크 완료** (위 Pre-Flight Check)
- [ ] Planning 결과(stages/03)에서 실험 설계 및 기준 확인
- [ ] Step 0: 데이터 접근성 검증 완료
- [ ] Step 1: 기존 분석 재현 완료
- [ ] Step 2: 결과 재현 확인 (수치 기준 충족)
- [ ] Step 3: 스몰 서브셋 모듈별 구현-검증 (A→B→C→D→E→F)
- [ ] Step 4: 베이스라인 비교 + Ablation 실험 완료
- [ ] Step 5: 실제 규모 데이터 벤치마크 완료
- [ ] 성능 미달 시: 원인 분석 및 인사이트 도출
- [ ] outputs/analysis/run_XX/에 코드, 결과, 로그 저장
- [ ] stages/04_analysis.md의 Current State 업데이트 (누적 경험 반영)

## Key Questions
- 타겟 성능을 달성했는가?
- 어떤 피처/모델이 효과적이었는가?
- 어떤 시도가 실패했고, 왜 실패했는가?
- 다음 시도에서 무엇을 바꿔야 하는가?

## 성능 미달 시
1. 결과를 사용자에게 보고
2. 원인 분석 결과 공유
3. 이전 단계로 돌아갈지, 계속 진행할지 사용자에게 확인
4. 루프백 판단 기준: docs/09_loopback_protocol.md 참조

## 산출물 예시
- `outputs/analysis/run_01/` — 코드, 모델, 결과, 로그

## Output Saving Guidelines
- 모든 중간 산출물은 outputs/ 아래에 파일로 저장 (채팅에만 작성 금지)
- 파일명 규칙: {주제}_{날짜}.md
- outputs/{stage}/run_{NN}/ 디렉토리 아래에 저장
- 코드 스니펫, 설정 파일, 로그도 동일하게 outputs/에 저장

## MANIFEST 인덱스 관리
- `outputs/MANIFEST.md` (전체 인덱스) + 각 하위 디렉토리별 `MANIFEST.md`
- 필수 항목: 파일 타입 태그(code/report/spec/checkpoint/result/log), 1줄 설명, 핵심 수치 요약
- 새 파일 추가 시마다 해당 디렉토리의 MANIFEST.md 업데이트
