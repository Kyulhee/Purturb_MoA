# Nexus Science — Windows 배포판 (Closed Beta)

NexusBrain 기반 AI 연구 에이전트의 **Windows 전용 단일 실행파일 배포**.
소스 설치, Node.js/Bun 설치, 의존성 설정 없이 바로 실행할 수 있습니다.

> **Closed beta**: 이 저장소는 초대받은 사용자에게만 공개됩니다. API
> 키나 피드백을 공유하려면 초대한 사람에게 연락하세요.

---

## 설치 (3분)

### 1. 저장소 복제

```powershell
git clone https://github.com/GooTec/nexus-science-win.git
cd nexus-science-win
```

또는 ZIP 다운로드: 저장소 페이지 → "Code" → "Download ZIP" → 압축 해제.

이 저장소에는 **README, LICENSE, .env.example, CHANGELOG** 만 들어
있습니다. 실제 실행파일은 아래 단계에서 별도로 다운로드합니다.

### 1-1. Windows 실행파일 다운로드

바이너리는 저장소가 아니라 **GitHub Releases** 에서 관리됩니다
(파일 크기 때문에). 원하는 버전을 받으세요:

- **최신 버전**: https://github.com/GooTec/nexus-science-win/releases/latest
  에서 `nexus-science.exe` 클릭 → 다운로드
- **특정 버전**: Releases 탭에서 해당 태그 선택
- **CLI 사용자**: `gh release download --repo GooTec/nexus-science-win v0.1.3 --pattern nexus-science.exe`

다운로드한 `nexus-science.exe` 를 저장소 폴더의 `bin\` 아래에
넣으세요:

```powershell
mkdir bin
move ~\Downloads\nexus-science.exe bin\nexus-science.exe
```

(경로는 예시. `bin\nexus-science.exe` 위치에 있으면 됩니다.)

### 2. `.env` 설정

`.env.example` 을 `.env` 로 복사한 뒤 내용을 편집합니다:

```powershell
copy .env.example .env
notepad .env
```

최소한 다음 값은 반드시 채워야 합니다:

| 키 | 설명 | 예 |
|---|---|---|
| `OPENAI_API_KEY` | vLLM 서버 또는 OpenAI 호환 엔드포인트 API 키 | `338af479...` |
| `OPENAI_BASE_URL` | 엔드포인트 URL (`/v1` 포함) | `https://example.trycloudflare.com/v1` |
| `OPENAI_MODEL` | 모델 ID | `glm-5.1-fp8` |

`.env.example` 은 GLM 5.1 기본 프로필이 이미 채워져 있으므로 보통 API
키만 바꾸면 됩니다.

### 3. 실행

```powershell
bin\nexus-science.exe --bare
```

정상 부팅되면 배너에 다음이 표시됩니다:

```
 Provider  OpenAI (custom)
 Model     glm-5.1-fp8
 Endpoint  https://...
```

프롬프트에 `/help` 를 치면 사용 가능한 슬래시 커맨드 목록이 나옵니다.

---

## `--bare` 플래그는 뭐고 왜 필요한가요?

현재 GLM 5.1 서버의 컨텍스트 한도가 **32,768 토큰**입니다. Nexus
Science 의 기본 시스템 프롬프트(도구 설명 + 워크플로 가이드 포함)는
이미 28,000 토큰에 가까워서, 대화 턴마다 입력 공간이 거의 남지
않습니다.

`--bare` 플래그는 "최소 시스템 프롬프트" 모드입니다:

- 자동 메모리 / 스킬 자동 로드 / 플러그인 동기화 / hook / LSP 등을
  건너뜁니다
- 시스템 프롬프트가 약 5,000 토큰으로 줄어듭니다
- 그 외 기능(도구 호출, 샌드박스, Task Gate 등)은 전부 그대로 작동

**서버가 128k 컨텍스트로 업그레이드된 이후에는** `.env` 에서
`CLAUDE_CODE_MAX_OUTPUT_TOKENS=4096` 줄과 `--bare` 플래그 둘 다
제거하고 실행하면 기본 기능이 전부 돌아옵니다. 이 저장소의 다음
릴리스 노트에서 알려드립니다.

---

## 사용 예시

```
> 해마 신경재생 관련 최신 논문 3편 찾아서 요약해줘

[agent calls PaperSearch with max_results=3]
[agent receives tool results]
[agent composes response with verified DOIs]
```

에이전트가 할 수 있는 일:

- PubMed / arXiv / Semantic Scholar 논문 검색
- 논문 PDF 다운로드 + 파싱 (PDF 저장소가 허용하는 경우)
- NexusBrain 지식베이스 조회 (환경 설정 필요)
- 가설 생성 → 실험 설계 → 결과 분석
- 로컬 파일 읽기 / 쓰기 / 편집 (Edit / Write / Read 도구)
- 셸 명령 실행 (Bash 도구 — 실행 전 항상 확인)
- 작업 분해 + Task Gate 추적

**할 수 없는 일** (혹은 의도적으로 막은 일):

- 인터넷 전체 쓰기 작업 (git push, 클라우드 배포 등 — 사용자가
  명시적으로 허용해야 함)
- 출처 없는 사실 날조 (과학 무결성 규범 탑재. 도구가 확인하지 못한
  정보는 "unverified" 레이블이 붙거나 답변 거부)

---

## 파일 구조

```
nexus-science-win/
├── README.md              — 이 파일
├── .env.example           — 환경 설정 템플릿
├── .env                   — (사용자가 만듦, git에 포함되지 않음)
├── CHANGELOG.md           — 릴리스 노트
├── LICENSE                — 라이선스 정보
└── bin/
    └── nexus-science.exe  — 단일 실행파일 (약 131 MB)
```

`.exe` 는 Bun 런타임을 내장한 자족 바이너리입니다. 추가 DLL 이나
node_modules 가 필요하지 않습니다.

---

## 업데이트

```powershell
git pull
```

또는 새 ZIP 을 받아서 기존 폴더의 `bin\nexus-science.exe` 만 덮어쓰기
해도 됩니다. `.env` 는 건드리지 마세요 (귀하의 설정이 날아갑니다).

---

## 문제 해결

### "This app can't run on your PC" / SmartScreen 경고

처음 실행 시 Windows Defender SmartScreen 이 "알 수 없는 게시자" 경고를
띄울 수 있습니다. 이는 코드 서명 인증서가 없는 단일 배포라서 그렇습니다
(closed beta 단계에서는 정상). **"추가 정보" → "실행"** 을 눌러
허용하세요.

### 배너에 `Provider Anthropic` 이 뜹니다

`.env` 에 `CLAUDE_CODE_USE_OPENAI=1` 이 누락됐습니다. 이 플래그가 없으면
`OPENAI_*` 값들이 무시되고 기본 Anthropic 경로로 갑니다. 추가한 뒤 다시
실행하세요.

### `400 Bad Request: maximum context length` 에러

`CLAUDE_CODE_MAX_OUTPUT_TOKENS=4096` 이 `.env` 에 있는지 확인하고,
실행 시 `--bare` 플래그를 붙였는지 확인하세요. 이 두 가지 없이 GLM 5.1
(32k 컨텍스트) 로는 거의 실행되지 않습니다.

### `.env` 가 안 읽힙니다

- `.env` 와 `bin\nexus-science.exe` 가 **서로 상위-하위 관계**여야
  합니다. 즉 `nexus-science-win\.env` + `nexus-science-win\bin\nexus-science.exe`
  이 정상. 다른 폴더에서 실행하면 다른 `.env` 를 찾을 수 있습니다.
- 줄 끝에 백슬래시 (`\`) 를 붙이지 마세요. `.env` 는 셸 문법이 아니라
  단순 key=value 포맷입니다.
- 값을 따옴표로 감싸지 마세요 (필요할 때만).

### 네트워크 타임아웃 / rate limit

- 사내 VPN 을 통해서만 접근 가능한 엔드포인트라면 VPN 상태를 확인.
- API 키가 만료되었거나 할당량이 소진됐을 수 있습니다. 초대한 관리자에게
  문의.

### 크래시 / 예상치 못한 종료

crash log 를 찾아서 리포트해주세요:

```powershell
bin\nexus-science.exe --bare 2> crash.log
```

그 다음 `crash.log` 와 재현 단계를 GitHub Issues 에 올려주세요.

---

## 피드백 / 이슈 리포트

이 저장소의 **Issues** 탭을 사용하세요. 초대된 사용자는 모두 issue 를
열 수 있습니다. 리포트할 때 포함하면 좋은 것:

1. 이 저장소의 릴리스 버전 (`CHANGELOG.md` 또는 `nexus-science.exe --version`)
2. 재현 단계
3. 기대한 결과 vs 실제 결과
4. 스크린샷 또는 터미널 출력

---

## 보안 / 프라이버시

- `.env` 의 API 키는 **로컬에서만** 사용됩니다. 바이너리에 하드코딩되지
  않습니다. 이 저장소를 포크하거나 다른 사람과 공유해도 귀하의 키는 함께
  유출되지 않습니다 (단, `.env` 자체를 commit 하지는 마세요).
- `bin\nexus-science.exe` 는 Nexus Science 소스의 일부 함수가 실행 중에
  콜아웃하는 엔드포인트(PubMed, arXiv, Semantic Scholar, 설정한 모델
  프로바이더)를 제외하면 외부와 통신하지 않습니다. 텔레메트리 / 분석 /
  사용자 추적은 비활성화돼 있습니다 (`no-telemetry-plugin` 으로 빌드).
- 소스 저장소는 내부 GitLab 에 있고, 본 배포 저장소는 빌드 결과물만
  포함합니다.

---

## 라이선스

Closed beta 단계. 별도 고지 전까지는 초대받은 사용자 본인 테스트
목적으로만 사용 가능. 재배포 / 역공학 / 제3자 제공 금지. 자세한 사항은
`LICENSE` 참조.
