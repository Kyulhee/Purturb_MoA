# cometspy 설치 및 dFBA 데모 실행 보고서

**프로젝트**: OCT LLM XAI - Phase 4 (dFBA 동적 시뮬레이션) 설계 참고
**날짜**: 2026-04-26
**작성**: NexusScience Agent

---

## 1. 설치 결과

### cometspy v0.6.3 설치 성공

```
pip install cometspy
```

**설치된 패키지**:
| 패키지 | 버전 | 비고 |
|--------|------|------|
| cometspy | 0.6.3 | Python COMETS 인터페이스 |
| cobra | 0.31.1 | 제약기반 모델링 (의존성) |
| optlang | 1.9.0 | LP 최적화 언어 |
| swiglpk | 5.0.13 | GLPK Solver 바인딩 |
| python-libsbml | 5.21.1 | SBML 파싱 |
| pandas | 2.3.3 | 이미 설치됨 |
| numpy | 2.4.4 | 이미 설치됨 |

**설치 확인**:
- `import cometspy`: 성공
- `__version__`: 모듈에 직접 속성 없음, `cometspy.comets.__version__` = "0.6.3"
- 모든 서브모듈 (`model`, `layout`, `params`, `comets`): 임포트 성공

---

## 2. COMETS Java 코어 설치 필요 여부

### 결론: **필수 설치 필요**

cometspy는 Python 래퍼(wrapper)이며, 실제 시뮬레이션 계산은 Java로 작성된 COMETS 엔진이 수행합니다.

**현재 환경 상태**:
- `COMETS_HOME` 환경변수: **설정되지 않음**
- Java 런타임: **PATH에 없음**
- COMETS JAR 파일: pip 패키지에 **포함되지 않음**

**COMETS Java 코어 설치 절차**:
1. Java JDK 11+ 설치
2. COMETS v2.12.5 다운로드 (https://github.com/segrelab/COMETS)
3. `COMETS_HOME` 환경변수를 COMETS 설치 디렉토리로 설정
4. (선택) Gurobi 설치 및 `GUROBI_HOME` 설정 (상용 Solver)
5. 또는 or-tools (무료 Solver, COMETS에 번들됨) 사용

---

## 3. dFBA 데모 실행 결과

### Python 측 기능 (Java 없이 작동)

| 기능 | 상태 | 비고 |
|------|------|------|
| cobra 모델 로드 | 성공 | `cobra.io.load_model("textbook")` 사용 |
| COMETS 모델 변환 | 성공 | E. coli core (95 반응, 72 대사물) |
| 2종 모델 생성 | 성공 | e_coli_core + ecoli_mutant |
| exchange 반응 열기 | 성공 | `model.open_exchanges()` |
| Layout 생성 | 성공 | 1x1 grid, 20개 대사물 |
| Media 설정 | 성공 | glucose, acetate, O2, NH4 등 |
| Parameters 설정 | 성공 | timeStep, maxCycles, Vmax, Km 등 |
| 모델 파일 쓰기 (.cmd) | 성공 | COMETS 포맷으로 파일 생성 가능 |

### 시뮬레이션 실행 (Java 필요)

```
KeyError: 'COMETS_HOME'
```

`cometspy.comets.comets.__init__()`에서 `os.environ['COMETS_HOME']`을 읽으려 시도하나,
환경변수가 설정되지 않아 `KeyError` 발생. **예상된 실패**.

---

## 4. 발생한 오류 및 한계점

### 오류 1: `cobra.test` 모듈 제거
- **원인**: cobra 0.31.1에서 `cobra.test` 모듈이 제거됨
- **해결**: `cobra.test.create_test_model("textbook")` 대신 `cobra.io.load_model("textbook")` 사용
- **영향**: cometspy 문서의 예제 코드가 구버전 cobra 기준으로 작성됨

### 오류 2: COMETS_HOME KeyError
- **원인**: cometspy.comets.comets.__init__()이 `os.environ['COMETS_HOME']`을 직접 접근
- **해결**: COMETS Java 코어 설치 후 환경변수 설정 필요
- **영향**: Java 코어 없이는 시뮬레이션 실행 불가

### 한계점 1: pandas 2.x 호환성 문제
- **원인**: `model.py`의 `add_signal()`, `add_multitoxin()` 메서드가 `DataFrame.append()` 사용
- **상태**: pandas 2.3.3에서 `DataFrame.append()`는 이미 제거됨
- **영향**: signaling 기능 사용 시 런타임 오류 발생 가능
- **회피**: 기본 dFBA에서는 signal 사용 안 함

### 한계점 2: Windows 지원
- Windows에서는 `comets_scr.bat` 스크립트를 통해 COMETS 실행 (comets.py 라인 396-400)
- classpath 테스트 루틴이 Windows에서는 스킵됨 (라인 271)
- 환경변수 설정이 Unix와 다를 수 있음

### 한계점 3: 상용 Solver 의존성
- 기본 optimizer가 GUROBI (상용 라이선스 필요)
- or-tools (무료) 대안 지원하나 설정 필요
- GLPK는 일부 기능 미지원

### 한계점 4: Python 3.13 호환성
- cometspy 자체는 Python 3.13에서 정상 임포트됨
- 그러나 cobra, pandas 2.x 등 종속성에서 deprecation 경고/오류 가능
- cometspy 테스트가 Python 3.13에서 공식 검증되지 않았을 가능성

---

## 5. Phase 4 설계를 위한 권고사항

### 옵션 A: cometspy + COMETS Java (공간 시뮬레이션 필요 시)
- 장점: 2D 공간 시뮬레이션, 생물막 모델링, 다종 상호작용 시각화
- 단점: Java 의존성, 설치 복잡성, Docker 권장
- 권장: Docker 컨테이너로 COMETS + Java 환경 구축

### 옵션 B: dfba-python (순수 Python dFBA)
- 장점: Java 불필요, cobra + scipy 기반, 설치 간편
- 단점: 공간 시뮬레이션 미지원, well-mixed 가정
- 권장: 비공간 dFBA에 적합

### 옵션 C: 커스텀 dFBA 루프
- cobra.Model.optimize() + Euler/RK 적분을 직접 구현
- 최대 유연성, 코드 제어 가능
- 검증된 시뮬레이션 엔진이 아닌 자체 구현에 따른 검증 부담

### 권장 경로
Phase 4에서 비공간(well-mixed) dFBA가 우선 필요한 경우 **옵션 B(dfba-python)** 로 시작,
향후 공간 시뮬레이션이 필요해지면 **옵션 A(cometspy + Docker)** 로 전환하는 단계적 접근을 권장합니다.

---

## 6. 참고 파일

- 데모 스크립트: `outputs/framing/run_01/cometspy_dfba_demo.py`
- 본 보고서: `outputs/framing/run_01/cometspy_dfba_report_20260426.md`
