# 세션 갈무리 요약 (2026-04-26)

## 수행된 작업

### Phase 3-4 기술적 실현성 심층 분석
- GEM->GNN 입력 변환: 이종 그래프(HGTConv) 방식 권장
- FBA 정답 데이터 생성: COBRApy 병렬 실행, 10,000조합 ~10-30분 추정
- Surrogate Model: GNN+XGBoost 하이브리드 + Active Learning
- TOPSIS 가중치 민감도: Entropy weight + Pareto front로 완화
- FLYCOP: 저장소 삭제 확인, COMETS v2.12.4를 대체로 추천
- dFBA 수치 안정성: BDF/Radau 솔버 + adaptive time-stepping
- NSGA-II 접종 비율 최적화: pymoo 활용, 대리 모델로 90% 비용 감소

### 핵심 결론
- Phase 3 실현성: MEDIUM-HIGH
- Phase 4 실현성: MEDIUM (FLYCOP 불가, COMETS 대체 필요)
- Phase 3-4 양방향 의존: Active Learning 루프로 연결

### CLAUDE.md 업데이트
- Output Saving Guidelines 추가: 모든 중간 산출물은 outputs/에 파일로 저장
- 채팅에만 작성하고 파일 미저장 금지 규칙 명시

## 산출물 위치
- 분석 보고서: outputs/planning/run_02/phase3_4_feasibility_analysis.md
- CLAUDE.md: 프로젝트 루트 (Output Saving Guidelines 추가됨)

## 다음 세션에서 재개할 작업
1. FLYCOP 원논문(Perez et al., 2018) 확보 및 기능 분석
2. COMETS/cometspy 실제 설치 및 dFBA 데모 실행
3. GNN 대사 네트워크 변환 프로토타입 구현
4. COBRApy 병렬 FBA 벤치마크 수행
5. Phase 3-4 통합 아키텍처 상세 설계

## 확인된 오픈소스 도구
- PyTorch Geometric: https://github.com/pyg-team/pytorch_geometric
- COBRApy v0.31.1: https://github.com/opencobra/cobrapy
- COMETS v2.12.4: https://github.com/segrelab/COMETS
- cometspy v0.6.1: https://github.com/segrelab/cometspy
- pymoo: https://github.com/anyoptimization/pymoo
- DeepChem: https://github.com/deepchem/deepchem
