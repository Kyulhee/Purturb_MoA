# Planning — Guide

## Pre-Flight Resource Check (권장)
계산 집약적 설계(FBA 병렬, GNN 학습, dFBA 시뮬레이션)를 포함하는 계획 수립 전:

```
python -c "import psutil; m=psutil.virtual_memory(); d=psutil.disk_usage('C:/'); print(f'RAM: {m.available/(1024**3):.1f}/{m.total/(1024**3):.1f}GB ({m.percent}%) | Disk: {d.free/(1024**3):.1f}/{d.total/(1024**3):.1f}GB ({d.percent}%) | CPU: {psutil.cpu_count(logical=False)}P/{psutil.cpu_count(logical=True)}L')"
```

- 가용 RAM < 4GB 또는 Disk < 10GB 시, 대규모 실험 설계를 사용자에게 보고 후 조정
- Multi-Agent로 병렬 탐색 시 최대 3개 에이전트 제한 (RAM 40GB 기준)

## Checklist
- [ ] Framing 결과(stages/02)에서 베이스라인/타겟 수치 확인
- [ ] 실험 설계 작성 (모델, 피처, 하이퍼파라미터 범위)
- [ ] 데이터 전처리 파이프라인 정의
- [ ] 교차 검증/평가 전략 수립
- [ ] 리스크 식별 및 대안 준비
- [ ] **사용자 컨펌 획득** (타겟 성능 + 실험 설계 승인)
- [ ] outputs/planning/run_XX/에 산출물 저장
- [ ] stages/03_planning.md의 Current State 업데이트

## Key Questions
- 어떤 모델/방법론을 사용할 것인가? 왜?
- 데이터 분할 전략은? (train/val/test, 교차검증)
- 핵심 하이퍼파라미터와 탐색 범위는?
- 실패 시 대안은 무엇인가?

## 중요
- Planning → Analysis 전환 전 **반드시 사용자 승인** 필요
- stages/03에 `✅ confirmed by user [날짜]` 기록

## 산출물 예시
- `outputs/planning/run_01/experiment_plan.md` — 실험 설계서
