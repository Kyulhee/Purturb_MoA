# Git Commit & Push Policy

## 커밋 시점
- 각 Stage의 산출물이 outputs/에 저장된 직후 커밋
- stages/ 문서 업데이트 후 커밋
- CLAUDE.md 업데이트 후 커밋 (다른 변경사항과 함께 가능)
- 대규모 코드/분석 산출물은 의미 단위로 분리 커밋

## 커밋 메시지 규칙
- `feat:` 새 산출물/기능 추가
- `docs:` 문서 업데이트
- `fix:` 수정/버그픽스
- `chore:` 설정/구조 변경

## 푸시 시점
- Stage 완료 시 사용자에게 푸시를 제안
- **사용자가 명시적으로 요청/승인 시 푸시** (자동 푸시 금지)
- 커밋은 자동 진행, 푸시는 사용자 승인 후 진행이 기본 원칙

## README.md 업데이트
- 푸시 전 최신 stages/ 문서를 꼼꼼히 읽고 README.md를 최신 상태로 업데이트
- README.md는 stages/의 핵심 내용을 외부 독자도 이해할 수 있도록 요약한 공개 문서

## .gitignore 준수
바이너리(.pt, .pkl, .h5, .zip 등), .env, bin/, .claude/, .nexus/ 등 미공개 항목은 항상 제외
