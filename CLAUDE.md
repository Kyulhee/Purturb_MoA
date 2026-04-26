# CLAUDE.md — Project Orchestrator

## Current Stage
`03_planning`

## Stage Map
| Stage | Name | Guide | State File | Outputs |
|-------|------|-------|------------|---------|
| 01 | literature_review | docs/01_literature_review.md | stages/01_literature_review.md | outputs/literature_review/ |
| 02 | framing | docs/02_framing.md | stages/02_framing.md | outputs/framing/ |
| 03 | planning | docs/03_planning.md | stages/03_planning.md | outputs/planning/ |
| 04 | analysis | docs/04_analysis.md | stages/04_analysis.md | outputs/analysis/ |
| 05 | interpretation | docs/05_interpretation.md | stages/05_interpretation.md | outputs/interpretation/ |

## Rules
1. 각 단계 진입 시 해당 stages/ 문서를 먼저 읽고 워크플로우를 따를 것
2. 산출물은 outputs/ 이외에 생성하지 말 것
3. Planning → Analysis 전환 시 반드시 사용자 컨펌을 받을 것
4. 기존 run 결과를 덮어쓰지 말 것 (run_01, run_02... 보존)
5. 루프백 시 stages/ 파일에 사유를 기록할 것

## stages/ 문서 관리 원칙
- **stages/ = 압축된 최신 지식**: 사용자가 개별 run을 읽지 않아도 핵심 내용 파악 가능해야 함
- **검증된 인사이트만 유지**: 모든 run에서 검증된 핵심 지식을 "검증된 핵심 지식" 섹션에 통합
- **과거 실패는 1줄 요약**: "Run 이력" 섹션에 최소한으로 기록, 세부 내용은 outputs/ 참조
- **매 업데이트 시 정제**: 새 run 결과를 반영할 때 불필요한 중복/과거 임시 기록은 제거
- **rules와 reports 분리**: stages/는 규칙과 검증된 지식, outputs/는 상세 보고서와 코드

## Output Saving Guidelines
- 모든 중간 산출물(분석 보고서, 설계 문서, 실험 결과)은 outputs/ 아래에 파일로 저장할 것
- 채팅에만 작성하고 파일로 저장하지 않는 것은 금지
- 파일명 규칙: {주제}_{날짜}.md (예: phase3_4_feasibility_20260426.md)
- outputs/{stage}/run_{NN}/ 디렉토리 아래에 저장
- 코드 스니펫, 설정 파일, 로그 등도 동일하게 outputs/에 저장
- CLAUDE.md는 프로젝트 공통 지침이므로 중간 산출물 보관 위치로 사용하지 말 것

## Resource Management Rules
1. **분석/데이터 처리 전 시스템 리소스 필수 체크**: RAM, Disk, CPU 상태를 먼저 확인 후 작업 시작
2. **RAM 가용 임계값**: 가용 RAM < 4GB 시 메모리 집약적 작업(FBA 병렬, GNN 학습) 금지 → 사용자에게 보고
3. **Disk 가용 임계값**: 가용 Disk < 10GB 시 대량 데이터 생성 작업 금지 → 사용자에게 보고
4. **Multi-Agent 동시 실행 제한**:
   - 최대 동시 백그라운드 에이전트: 3개 (RAM 40GB 기준)
   - RAM 32GB 이하 환경: 최대 2개
   - 각 에이전트 예상 메모리: 2-4GB (GNN/FBA), 0.5-1GB (문헌 검색)
5. **리소스 체크 명령**: 작업 시작 전 반드시 실행
   ```python
   python -c "import psutil; m=psutil.virtual_memory(); d=psutil.disk_usage('C:/'); print(f'RAM: {m.available/(1024**3):.1f}/{m.total/(1024**3):.1f}GB ({m.percent}%) | Disk: {d.free/(1024**3):.1f}/{d.total/(1024**3):.1f}GB ({d.percent}%) | CPU: {psutil.cpu_count(logical=False)}P/{psutil.cpu_count(logical=True)}L')"
   ```
6. **OOM 방지**: 대규모 데이터셋 로드 시 chunk 단위 처리, 전체 로드 금지
7. **임계값 초과 시**: 작업 중단, 사용자에게 리소스 부족 알림, 대안(샘플링, 청킹, 에이전트 수 감소) 제안

## Git Commit & Push Policy
1. **커밋 시점**:
   - 각 Stage의 산출물이 outputs/에 저장된 직후 커밋
   - stages/ 문서 업데이트 후 커밋
   - CLAUDE.md 업데이트 후 커밋 (다른 변경사항과 함께 가능)
   - 대규모 코드/분석 산출물은 의미 단위로 분리 커밋
2. **커밋 메시지 규칙**:
   - `feat:` 새 산출물/기능 추가
   - `docs:` 문서 업데이트
   - `fix:` 수정/버그픽스
   - `chore:` 설정/구조 변경
3. **푸시 시점**:
   - **사용자가 명시적으로 요청 시 푸시** (자동 푸시 금지)
   - 단, 사용자가 "커밋 및 푸시까지 진행해"라고 지시한 경우 즉시 푸시
   - 커밋은 자동 진행, 푸시는 사용자 승인 후 진행이 기본 원칙
4. **.gitignore 준수**: 바이너리(.pt, .pkl, .h5, .zip 등), .env, bin/, .claude/, .nexus/ 등 미공개 항목은 항상 제외

## How to Resume
1. 이 파일에서 current_stage 확인
2. 해당 stages/ 문서 읽기 → Current State 파악
3. docs/ 가이드 확인
4. outputs/ 기존 산출물 파악
5. 작업 재개
