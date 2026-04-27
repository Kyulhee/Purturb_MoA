# Literature Review Run 03: When Does GNN Embedding Provide Value Over Tabular Features for Metabolic Network Modeling?

**Date**: 2026-04-27
**Scope**: GNN vs tabular baselines for metabolic network prediction; surrogate models for FBA; transfer learning for metabolic models
**Context**: E2E pipeline where HGTConv produces 32-dim graph embeddings concatenated with 137-dim knockout mask for XGBoost growth rate prediction. XGBoost-only (R2=0.91) outperforms GNN+XGBoost (R2=0.82). Edge pretraining hurts.

---

## Executive Summary

The literature strongly supports our empirical finding: **when the knockout mask already fully encodes the prediction-relevant information, GNN embeddings are redundant and can hurt performance.** GNN value emerges specifically when (1) the graph topology encodes information absent from the input features, (2) the task requires generalization to unseen structures, or (3) transfer across models is needed. However, even in these favorable settings, GNN gains are modest (4-16% accuracy improvement over tabular baselines) and pretraining can be harmful if the pretext task does not align with the downstream objective.

---

## 1. GNN for Metabolic Network Prediction

### 1.1 GraphGDel: Constructing and Learning Graph Representations of Genome-Scale Metabolic Models for Growth-Coupled Gene Deletion Prediction

| Item | Detail |
|------|--------|
| **Authors** | Ziwei Yang, Takeyuki Tamura |
| **Year** | 2025 (revised 2026) |
| **arXiv** | 2504.06316 |
| **Code** | github.com/MetNetComp/GraphGDel |

**Key Contribution**: First systematic pipeline for constructing graph representations from constraint-based metabolic models (COBRA-style JSON with metabolites, reactions, gene-reaction rules). Combines graph topology with gene (amino acid) and metabolite (SMILES) sequence data via a multi-module deep learning framework.

**Architecture** (4 modules):
- **Meta-M**: LSTM autoencoder for metabolite (SMILES) representation
- **Gene-M**: LSTM autoencoder for gene (amino acid) representation
- **Graph-M**: GCN on metabolite graph (graph convolution + link prediction + graph encoding/decoding)
- **Pred-M**: Element-wise multiplication fusion + FC + sigmoid for binary gene deletion prediction; multi-task loss (classification + reconstruction + link prediction)

**Graph Representation**: Metabolite graph constructed from genome-scale metabolic model. Input: COBRA-style JSON with reactions containing id, metabolites, gene_reaction_rule. Pre-built for e_coli_core, iMM904, iML1515.

**GNN vs Tabular Baselines**:

| Baseline | e_coli_core | iMM904 | iML1515 |
|----------|-------------|--------|---------|
| Deep feedforward NN (no graph) | -14.04% | -16.26% | -13.18% |
| Sequence-learning baseline (no graph) | -6.17% | -4.96% | -5.31% |
| Topology-aware graph aggregation (same graph) | -5.10% | -4.36% | -4.70% |

**Key Insight**: The largest gains are over non-graph baselines (13-16%), but even the topology-aware graph baseline is consistently outperformed (4-5%). This shows **both the graph construction AND the learning framework contribute** -- simply having a graph representation is insufficient. Critically, the task here (predicting which genes to delete for growth-coupled production) inherently requires understanding gene-metabolite relationships that cannot be encoded in a simple feature vector, unlike our knockout-to-growth task where the mask fully specifies the perturbation.

**Relevance to Our Project**: HIGH. GraphGDel uses GCN on a metabolite graph from the same type of COBRA models we use. The key difference is their task (predicting gene deletion strategies) requires understanding metabolic topology, while our task (predicting growth from knockout mask) may not. Their finding that graph+sequence > sequence-only > flat-features suggests GNN value depends on whether the graph encodes task-critical information not in the features.

---

### 1.2 FluxGAT: Integrating Flux Sampling with Graph Neural Networks for Unbiased Gene Essentiality Classification

| Item | Detail |
|------|--------|
| **Authors** | Kieren Sharma, Lucia Marucci, Zahraa S. Abdallah |
| **Year** | 2024 |
| **arXiv** | 2403.18666 |
| **Venue** | NeurIPS 2024 Workshop on AI for New Drug Modalities |
| **Code** | github.com/kierensharma/FluxGAT |

**Key Contribution**: GNN model that predicts gene essentiality from graphical representations of flux sampling data, eliminating the need for predefined objective functions (which FBA requires and which introduce observer bias).

**Architecture**:
- Uses **Mass Flow Graph (MFG)** as graph representation
- Nodes = reactions; task = node-level classification (essential vs. non-essential)
- Graph constructed from flux sampling data (not FBA)
- Node features from reactant/product information
- Node labels from gene-protein-reaction (GPR) rules
- k-fold cross-validation

**GNN vs Baselines**: FluxGAT achieved **almost double the sensitivity of FBA** for gene essentiality prediction. However, no comparison with tabular ML baselines (Random Forest, XGBoost) was reported.

**Key Insight**: GNN provides value here specifically because the task (gene essentiality) cannot be formulated as a simple feature-prediction problem -- it requires understanding the flow of mass through the network. FBA's limitation is the need for a predefined objective, and flux sampling + GNN bypasses this. The graph structure IS the information source, not just supplementary.

**Relevance to Our Project**: HIGH. FluxGAT demonstrates a scenario where GNN adds value: when the task requires understanding network flow that tabular features cannot capture. Our knockout-to-growth task may not have this property because the knockout mask already determines the FBA solution. FluxGAT's approach of using flux sampling rather than FBA as input is interesting -- it suggests that if we encoded flux distributions as features, tabular methods might work even better.

---

### 1.3 GATTACA: Graph Neural Network-Based Reinforcement Learning for Controlling Biological Networks

| Item | Detail |
|------|--------|
| **Authors** | Andrzej Mizera, Jakub Zarzycki |
| **Year** | 2025 |
| **arXiv** | 2505.02712 |

**Key Contribution**: Deep reinforcement learning with GNN graph convolutions to control Boolean network models of gene regulatory/signaling pathway networks for cellular reprogramming. GNNs encode the biological system's structure and dynamics into latent representations for the DRL agent.

**Architecture**:
- Graph convolution operations as function approximator for DRL action-value function
- Novel control problem formulation for Boolean networks under asynchronous update mode
- Pseudo-attractor concept for scalability

**GNN vs Baselines**: Demonstrated "scalability and effectiveness" on large-scale biological networks. No direct comparison with tabular/non-graph RL agents reported.

**Key Insight**: GNN provides value because the control task requires understanding network structure (which genes regulate which) -- information that is inherently graph-structured. A tabular agent would need to learn these dependencies from data, while GNN encodes them as inductive bias.

**Relevance to Our Project**: MODERATE. Shows GNN value when the task is sequential decision-making over graph-structured state spaces, which is different from our one-shot prediction task.

---

## 2. Deep Learning Surrogate Models for FBA / Metabolic Flux Prediction

### 2.1 Integrative Genome-Scale Metabolic Modeling and Machine Learning for Biofuel Biomass Prediction

| Item | Detail |
|------|--------|
| **Authors** | Neha K. Nair, Aaron D'Souza |
| **Year** | 2026 |
| **Source** | arXiv search result |

**Key Contribution**: Combines Yeast9 genome-scale metabolic model with ML to predict biomass flux in S. cerevisiae. FBA generated 2,000 flux profiles under varying glucose, oxygen, ammonium conditions.

**Results**:
- **Random Forest**: R2 = 0.99989
- **XGBoost**: R2 = 0.999
- VAE identified 4 metabolic clusters
- SHAP revealed glycolysis, TCA cycle, lipid biosynthesis as key determinants
- Bayesian optimization yielded 12-fold biomass flux increase

**Key Insight**: **XGBoost and Random Forest achieve near-perfect R2 (>0.999) on FBA surrogate prediction** without any graph neural network component. This is the strongest evidence that when the input features (nutrient constraints) fully determine the FBA output, tabular methods are more than sufficient. No GNN was needed or would likely add value.

**Relevance to Our Project**: VERY HIGH. This is the closest analog to our problem: predicting FBA outputs from input features. The near-perfect tabular performance (R2=0.999) confirms that FBA surrogate prediction is fundamentally a tabular problem when inputs fully specify the constraint space. Our R2=0.91 with XGBoost-only is consistent with this finding (lower due to discrete knockout space vs continuous nutrient constraints).

---

### 2.2 Hybrid Physics-Informed Metabolic Cybergenetics: ML Surrogates Informed by FBA

| Item | Detail |
|------|--------|
| **Authors** | Sebastian Espinel-Rios, Jose L. Avalos |
| **Year** | 2024 |

**Key Contribution**: Hybrid physics-informed dynamic modeling framework connecting gene expression and cellular metabolism. Uses ML surrogates informed by FBA to embed metabolic network physics into simpler macro-kinetic process rates. Enables single-level optimization (vs bilevel) while preserving key metabolic network knowledge.

**Architecture**: ML surrogate replaces bilevel FBA optimization with single-level optimization. The surrogate maps enzyme levels to metabolic exchange fluxes, preserving the input-output relationship of FBA while being differentiable and faster to evaluate.

**Key Insight**: Rather than using GNN, they use a standard ML surrogate that replaces the FBA inner loop. The surrogate is informed by FBA physics but does not use graph structure. This approach works because the surrogate only needs to learn the FBA input-output mapping, not the graph topology.

**Relevance to Our Project**: HIGH. This is a practical demonstration that FBA surrogate models do not require graph structure -- the FBA input-output mapping can be learned by simpler ML models. The physics-informed approach (embedding FBA knowledge into the surrogate) is more effective than graph-based approaches for this type of problem.

---

## 3. Active Learning for Metabolic Engineering

### 3.1 ART: Automated Recommendation Tool for Synthetic Biology

| Item | Detail |
|------|--------|
| **Authors** | Tijana Radivojevic, Zak Costello, Kenneth Workman, Hector Garcia Martin |
| **Year** | 2019 |
| **Venue** | Metabolic Engineering (not arXiv) |

**Key Contribution**: ML and probabilistic modeling tool for systematic metabolic engineering. Uses sampling-based optimization to recommend strains for the next engineering cycle with probabilistic production-level predictions. Validated on biofuels, hoppy beer, and fatty acid production.

**Key Insight**: ART uses ensemble models (Random Forest, etc.) rather than GNNs. The recommendation loop is an active learning cycle: predict -> recommend -> experiment -> retrain. This is the standard approach in metabolic engineering -- tabular ML + uncertainty quantification + sequential design. No graph structure is needed because the design space (which genes to modify) is already tabular.

**Relevance to Our Project**: HIGH. ART represents the practical state-of-the-art for metabolic engineering active learning. It uses tabular features + ensemble models, consistent with our finding that XGBoost-only outperforms GNN+XGBoost. The active learning component (recommending next experiments) is what our Module C was designed to do, but ART's approach suggests tabular uncertainty estimates may be more reliable than GNN-embedding-based ones.

---

## 4. Heterogeneous Graph Neural Networks for Biological/Biomedical Networks

### 4.1 COMET: Comprehensive Metapath-based Heterogeneous Graph Transformer for Gene-Disease Association

| Item | Detail |
|------|--------|
| **Authors** | Wentao Cui, Shoubo Li, Chen Fang, et al. |
| **Year** | 2025 |
| **arXiv** | 2501.07970 |

**Key Contribution**: Metapath-based heterogeneous graph transformer for gene-disease association prediction. Uses 7 metapaths with transformer-based aggregation to capture global contexts and long-distance dependencies. BioGPT-initialized node features.

**Relevance**: MODERATE. Demonstrates heterogeneous graph transformer design for biomedical networks, but the task (gene-disease association from knowledge graphs) is fundamentally different from metabolic flux prediction. The metapath approach could inform how we encode different edge types in our metabolic graph.

---

### 4.2 HGTDR: Advancing Drug Repurposing with Heterogeneous Graph Transformers

| Item | Detail |
|------|--------|
| **Authors** | Ali Gharizadeh, Karim Abbasi, Amin Ghareyazi, et al. |
| **Year** | 2024 |
| **arXiv** | 2405.08031 |

**Key Contribution**: Three-step approach: (1) heterogeneous knowledge graph construction, (2) heterogeneous graph transformer network, (3) relationship scoring. Addresses information loss from converting heterogeneous data to homogeneous form.

**Key Insight**: Explicitly identifies the problem of converting heterogeneous biological graphs to homogeneous representations as losing information. However, the method only "performs comparably to previous methods" -- heterogeneous graph modeling does not automatically provide gains.

**Relevance**: MODERATE. Shows that heterogeneous graph transformers do not always outperform simpler approaches even on inherently heterogeneous biological data. The comparable performance suggests the graph structure may not always be the bottleneck.

---

### 4.3 Heterogeneous GNN for Breast Cancer Diagnosis from Histopathology

| Item | Detail |
|------|--------|
| **Authors** | Akhila Krishna K, Ravi Kant Gupta, Nikhil Cherian Kurian, et al. |
| **Year** | 2023 |
| **arXiv** | 2307.08132 |

**Key Contribution**: Heterogeneous GNN capturing spatial and hierarchical relations between cell and tissue graphs in histopathological images. Cross-attention vs transformer for relationship modeling.

**Key Insight**: Heterogeneous GNN achieves higher accuracy with fewer parameters than transformer-based state-of-the-art. GNN inductive bias is beneficial when the data has explicit hierarchical structure (cell -> tissue).

**Relevance**: LOW for metabolic networks, but demonstrates that heterogeneous GNNs can outperform transformers when the hierarchy is well-defined.

---

### 4.4 MOTGNN: Multi-Omics Integration with Tree-Generated Graph Neural Network

| Item | Detail |
|------|--------|
| **Authors** | Tiantian Yang, Zhiqian Chen |
| **Year** | 2025/2026 |
| **arXiv** | 2508.07465 |

**Key Contribution**: Three-component framework:
1. **XGBoost for omics-specific supervised graph construction** (using feature importance to define edges)
2. **Modality-specific GNNs for hierarchical representation learning**
3. **Deep feedforward network for cross-omics integration**

**Results**: Outperforms SOTA by 5-10% in accuracy, ROC-AUC, F1-score on three disease datasets. Robust to class imbalance. Built-in interpretability.

**Key Insight**: **XGBoost is used to CONSTRUCT the graph, not to replace it.** This is a hybrid approach where tabular methods define the graph structure and GNNs learn from it. The graph is not given a priori -- it is learned from data using XGBoost feature importance. This suggests a potential alternative for our project: use XGBoost feature importance to identify which gene-metabolite relationships matter for growth prediction, then construct a sparse GNN from these.

**Relevance**: HIGH. Demonstrates a principled way to combine XGBoost and GNN -- let tabular methods identify important relationships, then let GNN learn from the resulting structure. This hybrid approach could address our redundancy problem.

---

## 5. Transfer Learning and Pretraining for Graph Neural Networks

### 5.1 Toward a Universal Foundation Model for Graph-Structured Data

| Item | Detail |
|------|--------|
| **Authors** | Sakib Mostafa, Lei Xing, Md. Tauhidul Islam |
| **Year** | 2026 |
| **arXiv** | 2604.06391 |

**Key Contribution**: Graph foundation model using feature-agnostic structural prompts (degree statistics, centrality measures, community structure indicators, diffusion-based signatures) for cross-domain transfer learning. Pretrained once on heterogeneous graphs, then reused on unseen datasets with minimal adaptation.

**Results**: On SagePPI benchmark, supervised fine-tuning of the pretrained backbone reached **mean ROC-AUC of 95.5%**, a gain of **21.8% over the best supervised message-passing baseline**.

**Key Insights for When Pretraining Helps**:
1. **Feature-agnostic design is critical** -- relying on universal structural properties rather than node-specific features enables cross-domain transfer
2. **Heterogeneous graph pretraining transfers** -- pretraining on diverse graph types creates useful representations for unseen biomedical networks
3. **Large gains when supervised methods struggle** -- the 21.8% improvement suggests pretraining provides biggest advantage when supervised methods cannot learn generalizable representations from limited data
4. **Biomedical networks especially benefit** -- variability across cohorts/institutions makes transfer learning particularly valuable

**Relevance**: VERY HIGH. This paper directly addresses when GNN pretraining provides value: when the target domain has limited data AND when the pretraining uses feature-agnostic structural properties. Our edge prediction pretraining fails precisely because it uses node-specific features (metabolite/reaction/gene embeddings) rather than structural properties. A structural prompt approach could potentially work better.

---

### 5.2 DIB-OD: Preserving the Invariant Core for Robust Heterogeneous Graph Adaptation

| Item | Detail |
|------|--------|
| **Authors** | Yang Yan, Qiuyan Wang, Tianjin Huang, et al. |
| **Year** | 2026 |
| **arXiv** | 2604.10882 |

**Key Contribution**: Addresses negative transfer and catastrophic forgetting in GNN pretraining across heterogeneous domains. Decomposes representations into orthogonal invariant and redundant subspaces using Information Bottleneck distillation + HSIC. Self-adaptive semantic regularizer gates label influence based on predictive confidence.

**Key Insight on When Pretraining Hurts**:
- **Negative transfer**: occurs when methods cannot separate invariant (task-relevant) features from domain-specific noise, causing source-domain knowledge to actively harm target-domain performance
- **Catastrophic forgetting**: the invariant core gets overwritten during target-domain adaptation, especially when label influence is not properly controlled
- **Solution**: disentangle invariant vs. redundant subspaces; gate label influence based on confidence

**Relevance**: VERY HIGH. Our finding that edge prediction pretraining hurts (R2=0.91 vs 0.96) is a clear case of negative transfer. The edge prediction pretext task learns graph structure that is redundant with the knockout mask, creating domain-specific noise that harms the downstream growth prediction task. DIB-OD's framework explains WHY our pretraining fails and suggests a fix: disentangle structure learning (which may be useful for transfer) from feature learning (which may be task-specific and harmful to transfer).

---

## 6. GNN Limitations: When Simpler Methods Suffice

### 6.1 Machine Learning Small Molecule Properties in Drug Discovery (Review)

| Item | Detail |
|------|--------|
| **Authors** | Nikolai Schapin, Maciej Majewski, Alejandro Varela, Carlos Arroniz, Gianni De Fabritiis |
| **Year** | 2023 |
| **arXiv** | 2308.12354 |

**Key Finding**: "Neural networks, while more flexible, **do not always outperform simpler models**." Despite multiple diverse approaches, performances are "often comparable." The availability of high-quality training data remains crucial. Calls for standardized benchmarks and additional performance metrics.

**Relevance**: HIGH. This is a comprehensive review confirming that in molecular property prediction (the closest analog to metabolic network prediction), GNNs do not reliably beat simpler tabular methods. The key determinants of performance are data quality and quantity, not model complexity.

---

### 6.2 TrustworthyMS: Uncertainty-Aware Metabolic Stability Prediction with Dual-View Contrastive Learning

| Item | Detail |
|------|--------|
| **Authors** | Peijin Guo, Minghui Li, Hewen Pan, et al. |
| **Year** | 2025 |
| **Venue** | ECML-PKDD 2025 |

**Key Contribution**: Contrastive learning framework addressing two GNN limitations: (1) atom-centric message passing that disregards bond-level topological features, (2) lack of uncertainty quantification. Uses edge-induced feature propagation and Beta-Binomial uncertainty modeling.

**Key Insight**: Even within GNN design, standard message passing loses important information (bond-level features). Simply using a GNN is not enough -- the architecture must be carefully matched to the information structure of the problem.

**Relevance**: MODERATE. Reinforces that GNN architecture matters and standard message passing may not capture relevant information even when graph structure is beneficial.

---

### 6.3 GraphCliff: Short-Long Range Gating for Activity Cliffs

| Item | Detail |
|------|--------|
| **Authors** | Hajung Kim, Jueon Park, Junseok Choe, et al. |
| **Year** | 2025 |

**Key Finding**: Graph embeddings fail to separate similar molecules with large potency differences (activity cliffs). Standard GNNs over-smooth representations of structurally similar molecules, making them unable to capture subtle but critical differences.

**Relevance**: MODERATE. Analogous to our finding that GNN embeddings cannot distinguish between knockout masks that differ only in which genes are off -- the GNN smooths over these differences because the graph topology is the same regardless of which genes are knocked out.

---

## 7. Synthesis: When Does GNN Add Value Over Tabular Features?

### 7.1 Conditions Where GNN Adds Value

| Condition | Evidence | Mechanism |
|-----------|----------|-----------|
| **Graph topology encodes task-critical information not in features** | GraphGDel: +13-16% over flat features; FluxGAT: 2x sensitivity vs FBA | GNN captures relational structure (gene-metabolite, reaction-metabolite) that flat features miss |
| **Transfer to unseen organisms/models** | Graph Foundation Model: +21.8% over supervised baseline | Feature-agnostic structural prompts enable cross-domain transfer |
| **Task requires understanding network flow/dynamics** | FluxGAT, GATTACA | GNN encodes mass flow / regulatory dynamics as inductive bias |
| **Limited labeled data in target domain** | Graph Foundation Model, DIB-OD | Pretrained representations provide useful initialization when target data is scarce |
| **Graph structure must be discovered from data** | MOTGNN: XGBoost constructs graph, then GNN learns | Hybrid: tabular methods identify important relationships, GNN learns from them |

### 7.2 Conditions Where GNN Does NOT Add Value (or Hurts)

| Condition | Evidence | Mechanism |
|-----------|----------|-----------|
| **Input features already fully determine the output** | Nair & D'Souza: XGBoost R2=0.999 for FBA surrogate; Our result: XGBoost-only R2=0.91 | No information in graph topology that is not already in the feature vector |
| **Pretraining pretext task misaligned with downstream task** | Our edge pretraining: R2=0.91 vs 0.96; DIB-OD framework | Negative transfer: pretext task learns domain-specific noise rather than invariant structure |
| **Graph structure is identical across samples** | GraphCliff: GNN fails on activity cliffs; Our metabolic graph | GNN message passing produces similar embeddings for different knockout masks on the same graph |
| **Sufficient labeled data available** | Schapin et al. review: neural networks "do not always outperform simpler models" | With enough data, tabular methods learn feature-output mapping without needing structural inductive bias |
| **Discrete perturbations on a fixed graph** | Our knockout mask on textbook model | The graph does not change between samples; only the mask changes. GNN cannot differentiate samples based on graph structure. |

### 7.3 Key Insight for Our Project

The fundamental issue is that **our metabolic graph is static** -- the same graph (72 metabolites, 95 reactions, 137 genes) is shared across all samples. The only thing that varies between samples is which genes are knocked out (the mask). In this setting:

1. The GNN produces the same base embedding for every sample, with only minor variation from the knockout mask propagation through message passing
2. The knockout mask is already a complete description of the perturbation -- the FBA solution is fully determined by which genes are off
3. The GNN embedding is therefore redundant with the mask, and concatenating it introduces noise

**GNN would add value if**:
- The graph varied between organisms (different topology = different information)
- Some genes had unknown or partially known functions (GNN could infer from network context)
- The task required predicting properties that depend on pathway structure (e.g., "which metabolites accumulate?" rather than just growth rate)
- Transfer learning to a new organism where the knockout mask alone is insufficient

---

## 8. Recommendations

### 8.1 Short-term (Within Current Project)

1. **Accept XGBoost-only as the primary model** for growth rate prediction on the textbook model. The literature consistently shows that when input features fully determine FBA output, tabular methods are sufficient or superior.

2. **Abandon edge prediction pretraining**. The DIB-OD framework confirms that misaligned pretraining causes negative transfer. If pretraining is attempted, use feature-agnostic structural prompts (per the Graph Foundation Model approach).

3. **Reframe GNN value proposition**: Instead of "GNN improves prediction on a single model," test "GNN enables transfer across models." This is where the literature suggests GNN has the strongest case.

### 8.2 Medium-term (Design Modifications)

4. **MOTGNN-style hybrid**: Use XGBoost feature importance to construct a task-relevant graph (which gene-metabolite edges matter for growth?), then use GNN on this sparse graph. This addresses the redundancy problem by making the graph task-specific.

5. **Cross-model evaluation**: Test GNN+XGBoost on iJO1366 (1,367 genes) or iML1515. If GNN adds value on larger models, this would confirm the hypothesis that GNN helps when the feature space is too large for tabular methods alone.

6. **Structural prompts instead of learned embeddings**: Following the Graph Foundation Model approach, use fixed structural features (degree, centrality, community membership) rather than learned GNN embeddings. These are feature-agnostic and may transfer better.

### 8.3 Long-term (Research Direction)

7. **GNN for pathway-level predictions**: Instead of predicting a single scalar (growth rate), predict pathway fluxes or metabolite concentrations. These depend on network structure in ways that growth rate does not.

8. **Transfer learning across organisms**: Train GNN on multiple organisms' metabolic models and test generalization to unseen organisms. This is the strongest case for GNN value based on the literature.

9. **Active learning with structural uncertainty**: Instead of using GNN embedding-space uncertainty for AL, use structural uncertainty (e.g., which parts of the metabolic network are undersampled). This may be more reliable than embedding-based uncertainty.

---

## 9. Paper Summary Table

| # | Paper | Year | Graph Type | GNN > Tabular? | Conditions for GNN Value | Relevance |
|---|-------|------|-----------|----------------|-------------------------|-----------|
| 1 | GraphGDel (Yang & Tamura) | 2025 | Heterogeneous (metabolite graph) | YES (+13-16% over flat, +4-5% over graph baseline) | Graph topology + sequence data jointly needed | HIGH |
| 2 | FluxGAT (Sharma et al.) | 2024 | Homogeneous (MFG, reaction nodes) | YES (2x sensitivity vs FBA) | Task requires network flow understanding; FBA objective bias | HIGH |
| 3 | GATTACA (Mizera & Zarzycki) | 2025 | Homogeneous (Boolean network) | Not compared | Sequential control over graph-structured dynamics | MODERATE |
| 4 | Nair & D'Souza (Yeast9+ML) | 2026 | None (tabular) | N/A (tabular R2=0.999) | Tabular suffices when features determine output | VERY HIGH |
| 5 | Espinel-Rios & Avalos (Cybergenetics) | 2024 | None (ML surrogate) | N/A (no GNN) | FBA input-output mapping is tabular | HIGH |
| 6 | ART (Radivojevic et al.) | 2019 | None (ensemble ML) | N/A (tabular) | Active learning for strain design is tabular | HIGH |
| 7 | MOTGNN (Yang & Chen) | 2025 | Heterogeneous (XGBoost-constructed) | YES (+5-10% over SOTA) | XGBoost constructs graph; GNN learns from it | HIGH |
| 8 | Graph Foundation Model (Mostafa et al.) | 2026 | Heterogeneous (structural prompts) | YES (+21.8% over supervised) | Transfer learning with feature-agnostic prompts | VERY HIGH |
| 9 | DIB-OD (Yan et al.) | 2026 | Heterogeneous | YES (over SOTA transfer) | Disentangling invariant/redundant prevents negative transfer | VERY HIGH |
| 10 | Schapin et al. (Review) | 2023 | Various | NOT ALWAYS | Data quality matters more than model complexity | HIGH |
| 11 | COMET (Cui et al.) | 2025 | Heterogeneous (metapath) | YES (over SOTA) | Metapath-based aggregation for knowledge graphs | MODERATE |
| 12 | HGTDR (Gharizadeh et al.) | 2024 | Heterogeneous (KG) | COMPARABLE | Heterogeneous modeling does not automatically help | MODERATE |
| 13 | TrustworthyMS (Guo et al.) | 2025 | Homogeneous (molecular) | Not compared | GNN architecture must match problem structure | MODERATE |
| 14 | GraphCliff (Kim et al.) | 2025 | Homogeneous (molecular) | NO (fails on activity cliffs) | GNN over-smooths structurally similar samples | MODERATE |

---

## 10. Searched but Unavailable

The following repositories returned 404 and may be private, renamed, or removed:
- `github.com/mims-harvard/GemFL` -- GNN for metabolic flux (may exist under different org/name)
- `github.com/Chloroplast89/Meta-GNN` -- not found

No PyG examples specifically for metabolic/biochemical networks were found on the PyG documentation or tutorials.

---

*This review was compiled from arXiv search results, GitHub repository documentation, and existing project outputs. Web search (Google Scholar, Semantic Scholar) was unavailable during the review period. Some papers could not be fully verified against their PDFs due to arXiv ID mismatches or access limitations.*
