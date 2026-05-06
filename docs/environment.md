# Environment Setup — Guide

## 목적
Planning 단계 진입 시 프로젝트에 필요한 컴퓨팅 환경을 파악하고 기록한다.
실제 환경 정보는 `objects/current/environment.yaml`에 저장한다.

## Cold Start 대응 규칙

### AI 자체 설치 금지
AI가 임의로 Python, CUDA, 패키지를 설치하는 것을 금지한다. 이유:
- 사용자가 이미 구성한 환경(miniconda, venv 등)을 인지하지 못하고 중복/충돌 설치 가능
- GPU/CUDA 버전 불일치로 인한 런타임 오류 발생 위험
- 사용자 의도와 다른 CPU-only 패키지 설치 위험

### 환경 파악 순서
Planning 진입 시 아래 순서로 기존 환경을 먼저 파악:

1. **Python 런타임 탐지**: `where python`, `where python3`, `conda env list` 실행
2. **가상환경 확인**: conda env, venv, poetry 등 기존 환경 목록 파악
3. **CUDA/GPU 확인**: `nvidia-smi`, `python -c "import torch; print(torch.cuda.is_available())"` 로 GPU 가용성 확인
4. **핵심 패키지 버전**: 프로젝트에 필요한 주요 패키지의 설치 여부 및 버전 확인
5. **데이터 위치**: 프로젝트 데이터가 이미 로컬에 있는지, 다운로드가 필요한지 확인

### 파악 후 사용자와 확인
환경 파악이 완료되면 `objects/current/environment.yaml`에 초안을 작성하고 사용자에게 검증 요청:
- AI가 탐지한 환경이 실제 사용 의도와 일치하는지
- 어떤 런타임을 어떤 용도로 사용할지
- 누락된 환경이나 패키지가 있는지

## environment.yaml 구조

```yaml
runtime:
  - name: {환경명}
    path: {실행 경로}
    version: {Python 버전}
    usage: {사용 목적}
    cuda: {true/false}
    key_packages: [{주요 패키지 목록}]

invocation:
  {환경명}: {실행 명령어}

data_locations:
  {데이터명}: {경로}

hardware:
  gpu: {GPU 모델 및 VRAM}
```

## 업데이트 시점
- Planning 진입 시 (초안 작성)
- 새 런타임/패키지 설치 후
- 환경 문제로 인한 루프백 발생 시
