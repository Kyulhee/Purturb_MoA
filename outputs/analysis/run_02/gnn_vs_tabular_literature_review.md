# Literature Review: When Do GNNs Provide Value Over Tabular Features for Metabolic/Biochemical Network Modeling?

**Date:** 2026-04-27
**Context:** Experiment showed XGBoost-only (137-dim knockout mask) R2=0.91 vs GNN+XGBoost (32-dim GNN embedding + 137-dim mask) R2=0.82. GNN embedding was redundant because knockout mask already encodes gene perturbation info. Edge prediction pretraining was harmful (R2 dropped 0.96 to 0.91).

---

## Q1: When Does GNN Embedding Actually Add Predictive Value Over Raw Tabular Features in Biological/Metabolic Networks?

### Core Finding: GNNs Help When Graph Structure Encodes Information NOT Already Present in Node Features

The literature reveals a consistent pattern: **GNNs add value precisely when the prediction target depends on relational/positional information within the graph that is not captured by node-level features alone.** When node features already fully encode the relevant information (as in our knockout mask case), GNN message passing is redundant or harmful.

#### Conditions Where GNNs Add Value:

**1. Homophily-driven prediction tasks**
- Zhao et al. (2026, IEEE TNNLS): "Node Classification in GNNs: Impact of Neighborhood Label Distribution on Homophily and Heterophily" -- GNNs excel when connected nodes share labels/properties (homophily), because message passing effectively smooths predictions over the graph. When this assumption is violated (heterophily), GNNs can underperform feature-only MLPs.
- Chen et al. (2024, IEEE TNNLS 35(10)): "Exploiting Neighbor Effect: Conv-Agnostic GNN Framework for Graphs With Heterophily" -- Standard GNNs "perform well on homophilic graphs but may fail" on heterophilic ones due to the homophily assumption in graph convolution.
- **Implication for our case:** Metabolic network graphs are NOT homophilic -- adjacent reactions (sharing a metabolite) do NOT necessarily have similar essentiality. This explains why GNN message passing hurts rather than helps.

**2. When target depends on network position/flow, not just local features**
- Hasibi, Michoel & Oyarzun (2024, NPJ Syst Biol Appl): **FlowGAT** -- GNN outperformed SVC, MLP, and RF for gene essentiality prediction BECAUSE the Mass Flow Graph encodes directional flux propagation information. Flow Profile Embeddings (FPE) that capture inflow/outflow patterns up to k hops were critical -- they outperformed Local Degree Profile (LDP) and Random Walk Embedding (RWE). The graph structure here encodes *how metabolic flux propagates through the network*, which is NOT in any node feature.
- **Key insight:** FlowGAT works because the MFG edges carry FBA-derived mass flow weights that encode mechanistic information. Without these physics-informed edge weights, the graph topology alone would be far less informative.

**3. When the task requires aggregating context from neighbors**
- Alghamdi et al. (2021, Genome Res 31(10)): **scFEA** -- GNN is used to model metabolic flux balance constraints across a factor graph. The GNN architecture enforces that "the information aggregation between adjacent variables is constrained by the balance of the influx and outflux of each intermediate metabolite." Here, the graph structure IS the constraint -- flux balance literally requires neighbor information.
- Yuan & Nault (2025, Toxicol Sci 206(1)): **Metabolic network GNN for toxicant perturbation identification** -- GNN (3 GraphConv layers) applied to mouse Reactome pathways identified metabolic perturbations that pathway enrichment analysis missed (only 12/101 traditionally enriched reactions overlapped with GNN's top reactions). The GNN captured *hidden perturbation propagation* across the network.

**4. When features are noisy or incomplete and graph structure provides regularization**
- GNN message passing acts as a form of Laplacian regularization, smoothing noisy features over the graph. This is beneficial when features are unreliable but graph structure is trustworthy.

#### Conditions Where GNNs Are Redundant (Our Case):

- **Node features already encode the perturbation**: Our knockout mask (137-dim) already specifies which genes are perturbed. The metabolic graph adjacency does not add information about the perturbation itself.
- **Graph topology is known and fixed**: When the graph is a single static metabolic network (not varying across samples), it provides no discriminative information between samples. All samples share the same graph.
- **Task is not homophilic**: Gene essentiality / knockout effects do not necessarily propagate smoothly across metabolic graph neighborhoods.
- **No physics-informed edge weights**: If edges carry no quantitative flux or flow information, they only encode topology (which metabolites connect which reactions), which is already implicit in the stoichiometric matrix used by FBA.

---

## Q2: GNN as Surrogate Model for FBA / Constraint-Based Metabolic Modeling

### Papers Found:

**1. Song et al. (2025, Scientific Reports 15(1):6042)** -- "Coupling flux balance analysis with reactive transport modeling through machine learning"
- Uses **ANNs (not GNNs)** as reduced-order models (surrogates) for FBA within reactive transport simulations
- ANN-based surrogates achieved "substantial reduction of computational time by several orders of magnitude"
- Standard feedforward networks, not graph-structured
- **Key point:** For FBA surrogate modeling, the input is typically the environmental condition / medium composition, and the output is flux/growth. Graph structure is NOT needed because the surrogate is learning a mapping from condition-space to flux-space.

**2. Hasibi, Michoel & Oyarzun (2024)** -- FlowGAT (detailed above)
- Closest paper to "GNN surrogate for FBA" but with a crucial difference: it does NOT replace FBA with a GNN. Instead, it uses FBA solutions to CONSTRUCT the graph (Mass Flow Graph), then uses GNN for gene essentiality classification.
- The GNN learns from the *structure* of FBA solutions, not as a replacement for FBA.
- **Key insight:** GNN is useful as a *post-FBA analysis tool*, not as an FBA surrogate.

**3. Alghamdi et al. (2021)** -- scFEA (detailed above)
- GNN replaces iterative FBA solving with a learned model that enforces flux balance as a loss function
- The GNN is NOT a pure surrogate -- it is a *constrained neural optimizer* that uses the metabolic graph structure as the architecture itself
- This is the most creative use of GNN for metabolic modeling: the graph IS the model

**4. Espinel-Rios & Avalos (2024)** -- "Hybrid physics-informed metabolic cybergenetics"
- Uses ML surrogates informed by FBA to embed metabolic network physics into simpler macro-kinetic models
- Replaces bilevel optimization with single-level
- **Not GNN-based** -- uses standard neural networks

**5. Yang & Tamura (2025, IEEE TCBB 22(5))** -- DeepGDel
- Deep learning framework for predicting gene deletion strategies for growth-coupled production
- Uses "deep learning algorithms to learn and integrate sequential gene and metabolite data representation"
- 14.69%, 22.52%, and 13.03% accuracy improvements across three metabolic models
- **Architecture not specified as GNN** -- appears to use sequential/deep learning

### Key Insight for FBA Surrogate Modeling:
**No paper was found that uses GNN as a direct FBA surrogate.** The standard approach is:
- Input: environmental conditions, gene knockout status (tabular)
- Output: growth rate, flux distribution
- Model: MLP, ANN, or XGBoost
- **No graph structure needed** because the stoichiometric constraints are implicitly learned from training data

The graph is only useful when the TASK requires reasoning about network structure (e.g., which reactions to delete, how flux propagates, which metabolites are bottlenecks).

---

## Q3: Pretraining Strategies for GNN on Biochemical Networks

### Foundational Paper:

**Hu et al. (2020, ICLR Spotlight)** -- "Strategies for Pre-training Graph Neural Networks"
- **Critical finding:** Naive pre-training (graph-level OR node-level only) yields "limited improvement and can even lead to **negative transfer** on many downstream tasks."
- The proposed dual-level strategy (node-level + graph-level pretraining) "avoids negative transfer and improves generalization significantly" with up to 9.4% ROC-AUC improvement.
- **Directly relevant to our experiment:** Edge prediction pretraining (a form of node-level self-supervised pretraining) is exactly the type of naive pretraining that can hurt. When the downstream task (knockout effect prediction) requires different information than what edge prediction encodes, negative transfer occurs.

### Contrastive / Self-Supervised Approaches for Molecular/Biochemical GNNs:

**1. You et al. (2020, NeurIPS)** -- "Graph Contrastive Learning with Augmentations" (GraphCL)
- Four graph augmentation strategies for contrastive learning
- Achieves "similar or better generalizability, transferability, and robustness"
- **Key for biochemical graphs:** Choice of augmentation matters. Random node dropping or edge perturbation may destroy biologically meaningful structure.

**2. Veličković et al. (2019, ICLR)** -- "Deep Graph Infomax" (DGI)
- Maximizes mutual information between patch representations and high-level graph summaries
- "Does not rely on random walk objectives" (unlike DeepWalk/node2vec)
- Competitive or exceeding supervised performance on node classification benchmarks

**3. Xie et al. (2023, Brief Bioinform)** -- CAFE-MPP: "Self-supervised learning with chemistry-aware fragmentation for effective molecular property prediction"
- Fragment-based molecular graph under contrastive learning
- Chemistry-aware fragmentation respects molecular structure

**4. Jiang et al. (2024, Brief Bioinform)** -- DGCL: "Dual-graph neural networks contrastive learning for molecular property prediction"
- Dual-GNN contrastive learning with mixed molecular fingerprints
- Two-stage: pretraining then downstream feature training

**5. To et al. (2025, JCIM)** -- KGG: "Knowledge-Guided Graph Self-Supervised Learning to Enhance Molecular Property Predictions"
- Addresses "data scarcity and constrained model generalization" 
- Knowledge-guided self-supervised strategies

**6. He et al. (2024, Brief Bioinform)** -- POSIT: "Prototype-based contrastive substructure identification for molecular property prediction"
- Self-supervised framework for adaptively identifying meaningful substructures from molecular graphs

### What Would Work for Metabolic Network GNN Pretraining:

Based on the literature, **domain-appropriate pretraining** is essential:

1. **Knockout-aware pretraining** -- Pretrain on synthetic knockout-FBA pairs. The pretext task should predict flux/growth changes from knockout patterns, which is the actual downstream task. This is supervised pretraining, not self-supervised.

2. **Flux prediction as pretext** -- If using the metabolic graph, predict FBA-computed flux values from partial observations (mask some reaction fluxes, predict them from neighbors). This respects the physics of the network.

3. **Multi-task learning across conditions** -- Train on FBA solutions across many environmental conditions simultaneously. Different carbon sources, oxygen conditions, etc. create different flux distributions on the same graph, providing natural augmentation.

4. **AVOID: Generic link prediction / edge prediction** -- This was harmful in our experiment and is confirmed by Hu et al. as a source of negative transfer. Edge structure in metabolic networks is fixed and known; predicting it provides no useful inductive bias for knockout effect prediction.

---

## Q4: GNN for Metabolic Pathway Prediction, Gene Function Prediction, or Metabolic Engineering Optimization

### Metabolic Pathway Prediction:

**1. Baranwal et al. (2020, Bioinformatics 36(8))** -- "A Deep Learning Architecture for Metabolic Pathway Prediction"
- Early deep learning approach to predicting which metabolic pathways a compound participates in
- Achieved 96.08% accuracy
- Uses molecular structure features, not metabolic network graph structure

**2. Liu et al. (2024, IEEE/ACM TCBB 21(1))** -- MSGNN: "A Novel Multi-Scale Graph Neural Network for Metabolic Pathway Prediction"
- Multi-scale GNN with subgraph encoder, feature encoder, and global feature processor
- Accuracy, precision, recall, F1: 98.17%, 94.18%, 94.43%, 94.30%
- **Key point:** The graph here is the MOLECULAR graph (atoms and bonds), not the metabolic network graph. GNN is used to extract molecular structure features for pathway classification.

**3. Du et al. (2022, Bioinformatics 38(Suppl 1))** -- MLGL-MP: "Multi-Label Graph Learning Framework Enhanced by Pathway Interdependence for Metabolic Pathway Prediction"
- Leverages interdependencies between pathways for multi-label classification
- Compound encoder + pathway encoder + multi-label component

**4. Moozhippurath & Natarajan (2025)** -- "Enhancing Metabolomics Pathway Prediction with Sequential Graph Convolutional Network"
- Three-layer sequential GCN with ReLU activations on KEGG dataset
- Accuracy 98.00%, precision 92.10%, recall 93.02%
- Again uses MOLECULAR graphs (from SMILES), not metabolic network topology

**5. Hu et al. (2025, J Cheminform 17:56)** -- MotifMol3D: "Learning Motif Features and Topological Structure of Molecules for Metabolic Pathway Prediction"
- Combines motif information, GAT, and 3D structural data
- Uses XGBoost as final classifier (GAT for feature extraction + XGBoost for classification)
- Best precision (82.86%), recall (79.62%), F1 (81.21%)
- **Key finding:** Motif information captures biologically meaningful substructure patterns per pathway

**Critical observation:** All metabolic pathway prediction papers use GNN on MOLECULAR graphs (compound structure), not on metabolic NETWORK graphs (reaction topology). This makes sense because the question "which pathway does this compound belong to?" depends on the compound's chemical structure, not on its position in the metabolic network.

### Gene Essentiality / Function Prediction:

**1. Hasibi et al. (2024)** -- FlowGAT (detailed in Q2)
- GNN for gene essentiality using FBA-derived Mass Flow Graphs
- Outperforms SVC, MLP, RF with statistical significance
- **The only paper that truly uses metabolic network graph structure for gene essentiality prediction**

**2. Zhang et al. (2023, Brief Bioinform)** -- iEssLnc: "quantitative estimation of lncRNA gene essentialities with meta-path-guided random walks on the lncRNA-protein interaction network"
- Random-walk based (not GNN), on protein interaction network (not metabolic network)

### Metabolic Engineering Optimization:

**1. DeepGDel (Yang & Tamura, 2025)** -- Deep learning for gene deletion strategy prediction
- Not GNN-based; uses sequential deep learning
- Addresses growth-coupled production in genome-scale metabolic models

**2. SARTRE (Soleymani Babadi et al., 2023)** -- Prediction of metabolite-protein interactions using shadow prices from FVA
- Uses random forest on constraint-based features (shadow prices)
- "Highly competitive against recent deep-learning approaches"
- **Key insight:** Mechanistic features from constraint-based modeling (shadow prices) are powerful tabular features that capture network-level information without requiring GNN

### Phenotype Prediction:

**1. Triantafyllidis & Aguas (2025, NPJ Syst Biol Appl 11:92)** -- "Causality-aware graph neural networks for functional stratification and phenotype prediction at scale"
- Integrates MILP (CARNIVAL) for network reconstruction + GATv2Conv for classification
- Applied to TP53 mutation classification in CCLE (1,630 cell lines) and TCGA (12,471 tumors)
- Uses reconstructed gene regulatory networks as graph input
- TCGA misclassification: 16.17%; CCLE: 31.63%
- **Key innovation:** Engineered node features (color state, community-weighted centrality, in/out degree, Mode of Regulation) + "spotlight mechanism" for genes of interest

---

## Q5: The "Graph Structure Is Already Known" Problem

### The Core Issue:

When the graph topology is fixed and known (as in a genome-scale metabolic model), and node features already encode the perturbation (as in our knockout mask), the graph provides no discriminative information between samples. All samples share the same topology.

### What the Literature Says:

**1. Wu et al. (2019, IEEE TNNLS)** -- "A Comprehensive Survey on Graph Neural Networks"
- GNNs capture "dependence of graphs via message passing between the nodes of graphs"
- **The value comes from message passing, which requires that neighbors carry different information than the node itself.** If all nodes in a neighborhood have the same information (or information already captured by features), message passing is redundant.

**2. Xu et al. (2019, ICLR)** -- "How Powerful are Graph Neural Networks?"
- GNNs based on neighborhood aggregation "cannot learn to distinguish certain simple graph structures"
- GIN (proposed architecture) is "provably the most expressive among the class of GNNs" and equivalent to Weisfeiler-Lehman graph isomorphism test
- **Key implication:** If the prediction target can be expressed as a function of individual node features alone (without requiring relational reasoning), then the most expressive GNN is equivalent to an MLP applied to each node independently.

**3. Hu et al. (2020, ICLR)** -- Negative transfer from naive pretraining
- Node-level-only or graph-level-only pretraining can hurt downstream tasks
- This is exactly what happened in our experiment: edge prediction pretraining encoded irrelevant structural information that conflicted with the knockout prediction task

**4. Kundu et al. (2024, Biotechnology Advances)** -- "Machine learning for the advancement of genome-scale metabolic modeling"
- Review concludes that ML and CBM "mainly occur independently, which limits the concatenation of biological knowledge"
- The integration challenge is precisely about incorporating mechanistic knowledge (stoichiometry, flux balance) into ML models, not about using graph topology

### How Others Handle This Problem:

**1. Make the graph sample-specific (FlowGAT approach)**
- Hasibi et al. (2024) construct condition-specific Mass Flow Graphs where edge weights vary by carbon source
- Different conditions produce different graphs even on the same metabolic network
- This makes the graph structure informative: it encodes FBA-computed flux distributions that differ across conditions
- **Applicable to our case:** If we construct MFGs for each knockout condition, the graphs WOULD differ, and GNN could extract useful information from the flux propagation patterns

**2. Use the graph as architecture, not as input (scFEA approach)**
- Alghamdi et al. (2021) use the metabolic factor graph as the GNN architecture itself
- The graph structure constrains information flow, enforcing flux balance
- The graph is not "input data" but "model structure"
- **Applicable to our case:** Use the metabolic network to define message passing routes, where the GNN learns to predict fluxes subject to balance constraints

**3. Encode mechanistic information in edges/weights**
- FlowGAT: edge weights = FBA-computed mass flow
- scFEA: factor graph encodes stoichiometric relationships
- Without quantitative edge weights, the graph only encodes topology (which is known)
- **Applicable to our case:** If edges carry FBA-derived flux information, the graph becomes informative

**4. Use constraint-based features instead of graph structure**
- SARTRE (Soleymani Babadi et al., 2023): shadow prices from FVA are powerful tabular features
- Random forest on shadow prices outperformed deep learning on protein+metabolite features
- **Key insight:** The "graph information" can be captured as engineered tabular features (shadow prices, flux variabilities) without requiring GNN at all

**5. Multi-scale / hierarchical approaches**
- MotifMol3D (Hu et al., 2025): combines motif-level and atom-level features
- MSGNN (Liu et al., 2024): subgraph encoder + feature encoder + global processor
- These work on molecular graphs where multiple scales carry different information

---

## Synthesis: Why GNN Failed in Our Experiment and What Would Fix It

### Diagnosis:
1. **Static graph topology** -- All knockout samples share the same metabolic network. The graph provides zero discriminative power between samples.
2. **Feature sufficiency** -- The 137-dim knockout mask already fully encodes which genes are perturbed. The graph adds no new information.
3. **Harmful pretraining** -- Edge prediction pretraining forces the GNN to learn structural patterns (which reactions are connected) that are irrelevant to knockout effect prediction, causing negative transfer.
4. **Low embedding dimensionality** -- 32-dim GNN embedding was insufficient to capture the complex nonlinear relationship between gene knockouts and growth phenotypes, even if the graph had been informative.

### What Would Make GNN Useful:

| Approach | Why It Would Help | Reference |
|----------|------------------|-----------|
| Condition-specific MFGs (FlowGAT-style) | Graphs differ per knockout, encoding flux redistribution | Hasibi et al. 2024 |
| Factor graph with balance constraints (scFEA-style) | Graph IS the model architecture, not just input | Alghamdi et al. 2021 |
| FBA-derived tabular features (shadow prices, flux variabilities) | Captures network-level info as tabular features | Soleymani Babadi et al. 2023 |
| Pretrain on flux prediction, not edge prediction | Relevant pretext task avoids negative transfer | Hu et al. 2020 |
| Multi-condition training | Different conditions create different effective graphs | FlowGAT cross-condition results |

### Bottom Line:
**GNN is the wrong tool for our current task (predicting knockout effects from a knockout mask on a fixed metabolic network).** The literature consistently shows GNN adds value when:
1. Graph structure varies across samples or conditions
2. Target depends on relational/positional information not in features
3. Graph edges carry mechanistic/quantitative information
4. The task is homophilic (neighbors share properties)

None of these conditions hold in our setup. The correct approach is either:
- **Stick with XGBoost on tabular features** (including FBA-derived engineered features like shadow prices)
- **Redesign the GNN input** to use condition-specific Mass Flow Graphs where edge weights reflect FBA-computed flux distributions per knockout condition
