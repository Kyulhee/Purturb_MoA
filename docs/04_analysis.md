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

## Checklist
- [ ] **시스템 리소스 체크 완료** (위 Pre-Flight Check)
- [ ] Planning 결과(stages/03)에서 실험 설계 및 기준 확인
- [ ] 데이터 전처리 실행
- [ ] 모델 학습 및 평가
- [ ] 결과를 Planning에서 정의한 타겟과 비교
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
4. 루프백 시 stages/ 파일의 Previous Runs에 기록

## 산출물 예시
- `outputs/analysis/run_01/` — 코드, 모델, 결과, 로그
