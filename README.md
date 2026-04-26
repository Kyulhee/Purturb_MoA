# Purturb_MoA — Perturb-seq MoA Prediction Pipeline

대조 학습 기반 약물 임베딩으로 sci-Plex 단일 세포 섭동 데이터에서 MoA(Mechanism of Action) 클러스터링 품질을 개선하는 연구 프로젝트입니다.

## 프로젝트 개요

- **연구 질문**: "대조 학습 기반 약물 임베딩이 sci-Plex 단일 세포 섭동 데이터에서 MoA 클러스터링 품질을 개선하는가?"
- **데이터**: sci-Plex3 (188 compounds x 3 cell lines x 4 dosages = 649,340 cells, 7,561 genes, 16 MoA classes)
- **베이스라인**: GPAR-style DNN -> PANACEA -> chemCPA (r2=0.68 DEGs)

## 프로젝트 구조

```
Purturb_MoA/
├── CLAUDE.md                    # 프로젝트 오케스트레이터
├── docs/                        # 단계별 가이드
│   ├── 01_literature_review.md
│   ├── 02_framing.md
│   ├── 03_planning.md
│   ├── 04_analysis.md
│   └── 05_interpretation.md
├── stages/                      # 압축된 최신 지식 (단계별)
│   ├── 01_literature_review.md
│   ├── 02_framing.md
│   ├── 03_planning.md
│   ├── 04_analysis.md
│   └── 05_interpretation.md
└── outputs/                     # 단계별 산출물
    ├── literature_review/
    ├── framing/
    ├── planning/
    ├── analysis/
    └── interpretation/
```

## Phase 구성

### Phase 1-2: Perturb-seq MoA 예측 (기존)
- 대조 학습 + MoA-aware 임베딩
- 평가: leave-compound-out (분류), leave-MoA-out (클러스터링: silhouette, ARI, NMI)

### Phase 3: GNN+XGBoost 대리 모델 스크리닝
- GEM -> 이종 그래프 변환 (PyTorch Geometric HGTConv)
- FBA 정답 데이터 생성 (COBRApy 병렬)
- GNN+XGBoost 하이브리드 대리 모델
- Active Learning 루프

### Phase 4: dFBA 동적 시뮬레이션
- COMETS/cometspy 기반 미생물 군집 시뮬레이션
- NSGA-II 다목적 최적화 (접종 비율, pymoo)
- TOPSIS/Pareto 의사결정 지원

## 핵심 도구

| 도구 | 용도 |
|------|------|
| COBRApy | FBA (Flux Balance Analysis) |
| PyTorch Geometric | 이종 GNN |
| COMETS/cometspy | dFBA 동적 시뮬레이션 |
| pymoo | NSGA-II 다목적 최적화 |
| XGBoost | 대리 모델 회귀 |

## 핵심 설계 교훈 (이전 실패에서 학습)

| 결함 | 교훈 |
|------|------|
| Loss 불균형 (50:1) | effective gradient 크기로 weight 설계 |
| Drug encoder 동결 | fine-tuning 허용, LR 차등 적용 |
| leave-MoA-out 분류 평가 | 클러스터링 품질(silhouette, ARI, NMI)로 평가 |

## 라이선스

MIT
