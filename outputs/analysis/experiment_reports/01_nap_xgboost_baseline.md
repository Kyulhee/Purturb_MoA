# Experiment Report 01: NAP XGBoost Baseline -- GNN Embdedding Redundancy

**Run:** run_01
**Date:** 2026-04-26 ~ 2026-04-27
**Verdict:** FAIL (novelty insufficient)

---

## Hypothesis

GNN(HGT) embedding of metabolic network structure provides predictive value beyond raw tabular features (knockout mask) for FBA growth rate prediction. A GNN+XGBoost hybrid will outperform XGBoost-only.

**Origin:** Cross-domain novelty scan item -- "GNN surrogate for metabolic modeling" combining heterogeneous graph transformers with FBA surrogate modeling.

---

## Experimental Design

### Data
- Model: E. coli textbook (iML1515 subset), 137 genes, 72 metabolites, 95 reactions
- Heterogeneous graph: 3 node types (metabolite, reaction, gene), bidirectional edges (1,036 total)
- 837 FBA samples: 137 single-KO + 500 double-KO + 200 random (1-5 genes)

### Models Compared
| Model | Input | Dim | Notes |
|-------|-------|-----|-------|
| XGBoost-only | Knockout mask | 137 | Binary gene on/off vector |
| GNN(HGT)+XGBoost | GNN 32d embedding + mask | 169 | HGT encoder frozen after pretraining |
| No-pretrain GNN+XGBoost | GNN 153d embedding + mask | 290 | HGT encoder trained end-to-end |
| Edge-pretrain GNN+XGBoost | GNN 153d embedding + mask | 290 | HGT pretrained on edge prediction |

### Pretraining
- Edge prediction task: predict whether (metabolite, reaction) edge exists
- Transfer to downstream: freeze HGT encoder, concatenate embedding with mask for XGBoost

---

## Results

| Model | R2 (test) | Verdict |
|-------|-----------|---------|
| **XGBoost-only** | **0.9105** | Baseline |
| GNN+XGBoost (frozen, 32d) | 0.8236 | WORSE than baseline |
| No-pretrain GNN+XGBoost (153d) | 0.9596 | Slightly better (likely overfitting from higher dim) |
| Edge-pretrain GNN+XGBoost (153d) | 0.9058 | Worse than no-pretrain |

**Key numbers:**
- GNN embedding hurt: 0.9105 -> 0.8236 (delta = -0.0869)
- Pretraining hurt: 0.9596 -> 0.9058 (delta = -0.0538)
- Best GNN variant only marginally better than XGBoost-only, likely due to dimensionality increase (153d vs 137d)

---

## Failure Analysis

### Root Cause: Knockout Mask Already Encodes All Relevant Information

The 137-dimensional binary knockout mask is a **sufficient statistic** for FBA growth rate prediction. FBA is a deterministic function of which genes are knocked out -- the metabolic network topology is fixed and identical across all samples. Therefore:

1. **Graph topology is non-discriminative**: All samples share the same metabolic network. The graph provides zero discriminative power between samples.
2. **GNN message passing is noise**: Since the graph is static, HGT learns a fixed mapping from node features to embeddings. This mapping is redundant with the knockout mask itself.
3. **Edge prediction pretraining is misaligned**: Predicting edge existence (graph structure) is unrelated to predicting growth rate (perturbation effect). The pretrained weights encode structural information that interferes with downstream regression.

### Why This Was Not Obvious A Priori

- Literature shows GNNs add value when graph structure encodes relational information (FlowGAT, scFEA)
- Those cases have **varying graphs** or **physics-informed edge weights** (flux, flow)
- Our case has a **single static graph** with **no physics-informed edge weights** -- the exact conditions where GNNs are redundant

### Literature Confirmation (run_02)

The follow-up literature review (gnn_vs_tabular_literature_review.md) confirmed:
- GNNs help when prediction target depends on network position/flow NOT in node features
- GNNs are redundant when node features already encode perturbation (our case)
- Metabolic networks are NOT homophilic -- adjacent reactions don't share essentiality
- For FBA surrogate, input is condition-space, graph structure is not needed (Song et al., 2025)

---

## Knowledge Gained

1. **GNN value condition**: GNN embeddings are only useful when (a) graphs vary across samples, (b) edge weights encode physics/flux, or (c) node features are incomplete and graph structure provides regularization
2. **Pretraining alignment**: Pretraining task must be aligned with downstream task. Edge prediction != growth prediction
3. **XGBoost baseline is strong**: For tabular perturbation data with complete features, XGBoost is hard to beat
4. **Novelty gap**: "GNN surrogate for FBA" is not novel -- ANNs already serve this role (Song et al., 2025). The GNN angle adds nothing when the graph is static.

---

## Related Outputs

- `outputs/analysis/run_01/pipeline_e2e.py` -- E2E pipeline script
- `outputs/analysis/run_01/module_b_gnn_xgboost_surrogate.py` -- GNN+XGBoost module
- `outputs/analysis/run_01/pipeline_results.json` -- Raw results
- `outputs/analysis/run_01/phase3_4_technical_feasibility_report.md` -- Phase 3-4 analysis
- `outputs/analysis/run_02/gnn_vs_tabular_literature_review.md` -- Follow-up literature review
