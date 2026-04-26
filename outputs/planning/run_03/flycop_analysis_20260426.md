# FLYCOP 분석 보고서

> 작성일: 2026-04-26
> 원논문: Perez et al., BMC Systems Biology, 2018 (DOI: 10.1186/s12918-018-0639-6)
> 상태: 원논문 전문 확보 실패 (WebFetch/Springer/BioMed Central 접근 거부). 훈련 데이터 기반 분석이며, 실제 논문 확인 후 검증 필요.

---

## 1. FLYCOP 개요

**FLYCOP** (FuzzY Logic COmbined with Perturbation theory)는 미생물 군집의 계산적 설계를 위한 도구.
종 조합 및 초기 비율 최적화로 목표 생산물(바이오연료, 유기산 등) 극대화.

---

## 2. 핵심 알고리즘: Fuzzy Logic + Perturbation Theory 결합

### 2.1 Perturbation Theory 역할: "어디로 가야 하는지" 방향 제시
- 파라미터 섭동 후 시뮬레이션 응답으로 기울기 근사
- 현재 파라미터에서 각 차원별로 작은 변화(Δp)를 주고 결과 변화(Δf)를 관찰
- Δf/Δp로 국소 기울기 추정 → 다음 파라미터 이동 방향 결정

### 2.2 Fuzzy Logic 역할: "얼마나 좋은지" 평가
- 다목적 결과(생산량, 안정성, 자원 효율 등)를 membership function으로 통합
- 각 목적별로 "좋음/보통/나쁨" 정도를 0~1 fuzzy value로 표현
- Fuzzy aggregation으로 다목적 결과를 단일 스코어화

### 2.3 결합 방식
```
Perturbation 분석 → 기울기 근사 → 이동 방향 제시
       +
Fuzzy 평가 → 현재 상태의 "좋음" 정도 → 이동 크기 조절
       ↓
다음 파라미터 결정 → dFBA 실행 → 결과 반환 → 반복
```

---

## 3. 최적화 방법

- Grid/random search 대신 fuzzy-perturbation 결합으로 탐색 효율 향상
- 다목적(생산량 + 안정성)을 fuzzy aggregation으로 단일 스코어화 후 최적화
- 지역 탐색 위주 (전역 최적해 보장 불가)

---

## 4. dFBA 연동 방식

```
FLYCOP이 파라미터 설정 (종 비율, 환경 조건)
       ↓
dFBA가 GSMM 기반 동적 시뮬레이션 수행 (COBRApy 기반 추정)
       ↓
결과를 FLYCOP으로 반환 (생체량, 대사물 농도, 생산량)
       ↓
Fuzzy 평가 + Perturbation 분석 → 다음 파라미터 결정
       ↓
반복 (수렴 시 종료)
```

---

## 5. Python 3 호환성

- 2018년 발표 당시 Python 2/3 전환기
- COBRApy는 Python 3 지원 → 호환성 문제는 제한적
- 단, 구버전 라이브러리 의존성 확인 필요
- **GitHub 저장소 삭제됨** → 실제 코드 확인 불가

---

## 6. 주요 한계점

| 한계 | 상세 | 우리 프로젝트 영향 |
|------|------|-------------------|
| 1차 섭동 근사 | 비선형 영역에서 오차 큼 | Active Learning으로 대체 가능 |
| Fuzzy 규칙 주관성 | 멤버십 함수 설계자 의존 | TOPSIS/Entropy weight로 객관화 |
| 지역 최적해 함정 | 전역 최적해 보장 불가 | NSGA-II (다목적)로 완화 |
| well-mixed 가정 | 공간 구조 무시 | COMETS로 공간 시뮬레이션 가능 |
| 유전자 조절 누락 | 정적 GPR만 반영 | 시간 의존적 조절은 별도 모델링 필요 |
| 소규모 검증 | 2-3종만 검증 | 대규모 확장 시 추가 검증 필요 |
| Wet-lab 검증 제한적 | 실험 검증 결과 부족 | in silico 결과의 생물학적 타당성 주의 |

---

## 7. 우리 프로젝트와의 관련성

**낮음**: FLYCOP의 "perturbation theory"는 파라미터 공간 탐색 기법이며, 우리의 약물 섭동 예측(perturbation prediction)과는 개념이 다름.

**그러나 다음 측면에서 참고 가치 있음**:
1. 다목적 평가 + fuzzy 통합 → 우리의 TOPSIS 설계에 참고
2. dFBA + 최적화 루프 → 우리의 Active Learning 루프와 유사 구조
3. 종 비율 최적화 → Phase 4 NSGA-II 접종 비율 최적화와 동일 문제

---

## 8. 대체 전략 (FLYCOP 불가용)

| FLYCOP 기능 | 대체 도구 | 비고 |
|-------------|----------|------|
| dFBA 실행 | COMETS/cometspy | 공간 지원, 활발 유지보수 |
| 비공간 dFBA | dfba-python | 순수 Python, scipy 기반 |
| 파라미터 탐색 | pymoo (NSGA-II) | 다목적, 범용 |
| 다목적 평가 | TOPSIS + Entropy weight | 객관적 가중치 |
| Fuzzy 통합 | scikit-fuzzy | 필요시 사용 가능 |

---

## 9. 후속 작업

1. [ ] 웹 접근 복구 후 Springer/BioMed Central에서 원논문 전문 확인
2. [ ] 올바른 PMC ID 재검색 (PMC6201919는 타 논문)
3. [ ] FLYCOP 핵심 기능(다목적 fuzzy 평가)을 TOPSIS로 대체 구현
4. [ ] COMETS + NSGA-II 파이프라인으로 FLYCOP 기능 재현
