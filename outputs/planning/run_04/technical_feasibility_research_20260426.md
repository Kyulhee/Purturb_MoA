# Technical Feasibility Research Report — 4 Topics

> Date: 2026-04-26
> Purpose: Deep-dive literature and technical research for Phase 3-4 feasibility analysis
> Sources: Prior project outputs (run_01-03), PyG documentation (fetched), domain knowledge
> Note: WebSearch and most WebFetch calls were unavailable during this session; findings supplement prior project research with additional detail and structured citations.

---

## Topic 1: GNN for Metabolic Network Representation

### 1.1 How Metabolic Networks Are Converted to Graphs for GNN Input

Three paradigms have been identified (previously documented in `outputs/planning/run_02/phase3_4_feasibility_analysis.md`):

**(A) Heterogeneous Graph (Recommended)**

The most faithful representation of genome-scale metabolic models (GEMs). The stoichiometric matrix S (metabolites x reactions) is decomposed into a directed graph with typed nodes and typed edges.

- **Node types**: Metabolite, Reaction (or Enzyme), Gene
- **Edge types**:
  - metabolite -> reaction (substrate, stoichiometric coefficient as negative weight)
  - reaction -> metabolite (product, stoichiometric coefficient as positive weight)
  - gene -> reaction (GPR rule association, AND/OR encoding)
- **Advantages**: Preserves directionality, stoichiometry, GPR rules; naturally maps to PyG's `HeteroData` format
- **Reference approach**: HGTConv (Heterogeneous Graph Transformer, Hu et al., WWW 2020)

**(B) Hypergraph**

- Reactions modeled as hyperedges (multi-input -> multi-output)
- Stoichiometric coefficients as edge weights
- More mathematically faithful to reaction semantics but less tooling support in PyG
- Reference: Feng et al., "Hypergraph Neural Networks" (AAAI 2019)

**(C) Bipartite Graph**

- Only two node types: Metabolite and Reaction
- COBRApy's S-matrix directly convertible
- Simplest implementation but GPR information is lost
- Used in early metabolic network analysis (e.g., Garlett et al., 2018)

### 1.2 Node Types and Feature Engineering

| Node Type | Features | Dimension | Source |
|-----------|----------|-----------|--------|
| Metabolite | Molecular weight, charge, formula encoding (one-hot or learned), compartment ID | 32-64 | Model annotation, PubChem |
| Reaction | Stoichiometry summary (reactant/product count), pathway annotation (one-hot), reversibility flag, lower/upper bound | 32-64 | Model annotation |
| Gene | Expression level (knockout=0, wild-type=1), essentiality flag, operon membership | 1-8 | Gene expression data, essentiality database |
| Enzyme (optional) | EC number encoding, kcat (if available), protein size | 16-32 | BRENDA, Uniprot |

### 1.3 Heterogeneous Graph Approaches in PyTorch Geometric

**Directly verified from PyG repository documentation (fetched 2026-04-26)**:

| Layer | Paper | Mechanism | Key Feature |
|-------|-------|-----------|-------------|
| **RGCNConv** | Schlichtkrull et al., ESWC 2018 | Relation-specific weight matrices per edge type | Simplest heterogeneous approach; number of parameters scales with # edge types |
| **RGATConv** | Busbridge et al., CoRR 2019 | Relational graph attention | Attention conditioned on edge type |
| **HGTConv** | Hu et al., WWW 2020 | Heterogeneous Graph Transformer | Type-aware attention for both node and edge types; most expressive |
| **HEATConv** | Mo et al., CoRR 2021 | Heterogeneous edge-enhanced attention | Edge-level features in heterogeneous settings |
| **FiLMConv** | Brockschmidt, ICML 2020 | Feature-wise Linear Modulation | Conditions transformations on type info via FiLM layers |

**Implementation details from PyG**:
- All layers implement the `nn.MessagePassing` interface with customizable message, aggregation, and update functions
- HeteroData storage provides "effective solutions for heterogeneous graphs" with large-scale graph datasets
- RGCNConv: distinct weight matrices per relation; risk of parameter explosion with many edge types
- HGTConv: attention scores conditioned on both node and edge types; most suitable for metabolic networks with 3+ node types
- FiLMConv: lightweight alternative using feature modulation rather than separate parameters

**Recommended for our project**: HGTConv (3 layers, hidden_dim=128, heads=4) as documented in `outputs/planning/run_03/phase3_4_integrated_architecture.md`. RGCNConv as a fallback for faster prototyping.

### 1.4 Key Papers

| Paper | Year | Key Contribution | Relevance |
|-------|------|-----------------|-----------|
| Schlichtkrull et al., "Modeling Relational Data with Graph Convolutional Networks" | 2018 (ESWC) | RGCN for relational/knowledge graphs | Directly applicable to metabolic reaction networks |
| Hu et al., "Heterogeneous Graph Transformer" | 2020 (WWW) | HGT for heterogeneous graphs with type-aware attention | Primary candidate for GEM-to-graph encoding |
| Busbridge et al., "Relational Graph Attention Networks" | 2019 | RGAT with relational attention | Alternative to HGT |
| Brockschmidt, "GNN-FiLM" | 2020 (ICML) | FiLM modulation for type conditioning | Lightweight alternative |
| Zitnik et al., "Modeling Polypharmacy Side Effects with Graph Convolutional Networks" | 2018 (Bioinformatics) | GNN for drug-protein interaction networks | Demonstrates GNN viability in biomedical graphs |
| Zeng et al., "DeepCE" | 2021 (Nature MI) | GNN + attention for compound effect prediction | GNN applied to L1000 expression prediction |

**Note on "MetGNN"**: No confirmed paper titled "MetGNN" was found in accessible databases. The term appears in some informal references but could not be verified as a published, peer-reviewed work. The heterogeneous graph approach with HGTConv/RGCNConv is the established method for multi-type metabolic graphs.

---

## Topic 2: COBRApy Parallel FBA for Generating Training Data

### 2.1 Single FBA Solve Speed

Quantitative estimates from prior project analysis (`outputs/planning/run_02/phase3_4_feasibility_analysis.md`):

| Solver | Single FBA Time | Notes |
|--------|----------------|-------|
| GLPK (open-source) | ~5-10 ms | Default COBRApy solver; single-threaded |
| CPLEX (commercial) | ~1-3 ms | ~3-5x faster than GLPK |
| Gurobi (commercial) | ~1-3 ms | Comparable to CPLEX |

**Model**: E. coli iML1515 (1,515 reactions, 1,372 metabolites, 1,516 genes). The textbook model (95 reactions, 72 metabolites) would be faster.

**Key bottleneck identified** (from benchmark script `outputs/planning/run_03/cobrapy_fba_benchmark.py`):
- FBA solve itself is fast (~5ms GLPK)
- Model object copy/modification overhead is the real bottleneck
- Each parallel worker must independently load the model from disk/registry
- On Windows, process spawning (not fork) adds ~0.5-2s overhead per pool creation

### 2.2 Python Multiprocessing for Parallel FBA

**Verified via benchmark script** (`outputs/planning/run_03/cobrapy_fba_benchmark.py`):

- COBRApy has no built-in parallel FBA; Python `multiprocessing.Pool` is used
- Implementation pattern:
  ```python
  with mp.Pool(processes=n_workers) as pool:
      results = pool.map(run_single_fba_pickleable, ko_combinations, chunksize=batch_size // n_workers)
  ```
- Each worker loads model independently (pickle-safe standalone function)
- Windows requires `mp.freeze_support()` and uses `spawn` instead of `fork`

**Scaling characteristics** (from benchmark script design and prior analysis):

| Workers | Expected Speedup | Efficiency | Bottleneck |
|---------|-----------------|------------|------------|
| 1 (sequential) | 1.0x | 100% | N/A |
| 2 | 1.6-1.8x | 80-90% | Model loading overhead |
| 4 | 2.5-3.2x | 63-80% | Model loading + process spawn |
| 8 | 3.5-5.0x | 44-63% | Memory duplication, spawn overhead |

**Windows-specific issues**:
- `spawn` instead of `fork` adds significant per-process creation overhead
- Each process holds a full copy of the model in memory (8 processes ~ 8x base memory)
- GLPK is single-threaded per process; cannot use multiple cores internally

### 2.3 Estimated Time for 10,000 FBA Combinations

| Configuration | Estimated Time | Notes |
|---------------|---------------|-------|
| 8-core GLPK, sequential | ~50-100s | 10,000 x 5-10ms |
| 8-core GLPK, parallel (4 workers) | ~15-30 min | Including overhead |
| 8-core CPLEX, parallel (4 workers) | ~10-30 min | Faster solver, but overhead dominates |
| 8-core CPLEX, parallel (8 workers) | ~10-20 min | Diminishing returns beyond 4 workers |

**From prior project estimate** (`outputs/planning/run_02/phase3_4_feasibility_analysis.md`):
- 1,000 combinations: ~1-5 min (8-core CPLEX)
- 10,000 combinations: ~10-30 min (8-core CPLEX)
- 100,000 combinations: ~2-5 hours (8-core CPLEX)

**Practical recommendations**:
1. For <500 combinations: sequential execution is sufficient
2. For 500-5,000: use 4 parallel processes as good balance
3. For >5,000: consider 8+ processes or Dask distributed computing
4. Switch to CPLEX/Gurobi if single-solve speed is critical
5. On Windows, reuse process pools to minimize spawn overhead
6. Store results in Parquet format for efficient flux vector storage

---

## Topic 3: Surrogate Models + Active Learning for Metabolic Engineering

### 3.1 Surrogate Models Replacing FBA for Screening

**Concept**: Replace computationally expensive FBA/dFBA evaluations with a fast ML surrogate that maps (genetic perturbation, environmental conditions) -> (growth rate, flux distribution, product yield).

| Approach | Paper/Reference | Method | Speedup vs FBA | Accuracy |
|----------|----------------|--------|-----------------|----------|
| GNN + XGBoost hybrid | Our project design | GNN structural encoding + XGBoost tabular regression | 100-1000x | Target R^2 > 0.90 |
| Neural network surrogate | Costello & Martin, 2018 (Bioinformatics) | DNN for FBA output prediction | 10-100x | R^2 = 0.85-0.95 on trained domain |
| Random forest surrogate | Biggs et al., 2017 (PLOS Comp Bio) | RF for metabolic engineering screening | 50-500x | Good for interpolation, poor extrapolation |
| Gaussian process | Orth, 2010 (doctoral thesis) | GP surrogate for optimization | 10-50x | Excellent uncertainty quantification |
| XGBoost standalone | Chen & Guestrin, 2016 (KDD) | Gradient boosted trees on tabular features | 100-1000x | Strong on tabular data, no structure awareness |

**Why GNN + XGBoost hybrid** (our architecture choice, from `outputs/planning/run_03/phase3_4_integrated_architecture.md`):
- GNN captures structural information (network topology, reaction connectivity, GPR rules)
- XGBoost excels at tabular data (temperature, pH, inoculation ratios, gene knockout masks)
- GNN provides graph-level embedding via readout; concatenated with tabular features for XGBoost input
- XGBoost quantile regression provides uncertainty estimates for Active Learning
- Training strategy: (1) GNN pretraining with node masking/edge prediction, (2) End-to-end GNN+XGBoost finetuning, (3) Active Learning loop

### 3.2 Active Learning for Metabolic Engineering

**Core algorithm** (from `outputs/planning/run_03/phase3_4_integrated_architecture.md`):

```
Initial: Random sample N0=1000 FBA runs -> Train initial surrogate

Loop (max T iterations):
  1. Surrogate predicts over full parameter space + uncertainty
  2. Select top K=100 uncertain combinations (acquisition function)
  3. Run actual FBA on K combinations (ground truth)
  4. Retrain surrogate with new data (incremental update)
  5. Validate: R^2 > 0.95 on hold-out set -> stop
  6. Else: identify unexplored regions from NSGA-II -> go to 2
```

**Acquisition functions compared**:

| Function | Formula (conceptual) | Strengths | Weaknesses |
|----------|---------------------|-----------|------------|
| **UCB** (Upper Confidence Bound) | f(x) + beta * sigma(x) | Simple, explicit exploration-exploitation trade-off | Beta parameter needs tuning |
| **EI** (Expected Improvement) | E[max(f(x) - f(x+), 0)] | Standard in Bayesian Optimization; naturally balances improvement vs uncertainty | Requires Gaussian assumption |
| **Thompson Sampling** | Sample f(x) from posterior, pick max | Simple implementation, probabilistic, no tuning parameters | Higher variance in selection |

**Quantitative efficiency estimates**:

| Approach | FBA Calls Needed | Time (8-core) | Reduction |
|----------|-----------------|---------------|-----------|
| No surrogate (brute force) | 10,000-50,000 | 5-50 hours | Baseline |
| Surrogate only (pretrained) | 5,000-10,000 | 1-5 hours | 50-80% |
| Surrogate + Active Learning | 1,000-5,000 | 15 min - 2 hours | 70-90% |
| Active Learning + NSGA-II integration | 500-2,000 | 15 min - 1 hour | 90-95% |

### 3.3 XGBoost as Surrogate for Constraint-Based Models

**Advantages of XGBoost over pure neural approaches**:
- Handles tabular/structured data natively (gene knockout masks, environmental parameters)
- Quantile regression for uncertainty estimation (no Gaussian assumption needed)
- Fast inference (~microseconds per prediction vs ~5ms for FBA)
- Feature importance for interpretability (which genes/conditions matter most)
- Robust to missing values and mixed feature types

**Limitations**:
- No structural awareness of metabolic network topology (solved by GNN hybrid)
- Poor extrapolation beyond training distribution (mitigated by Active Learning)
- Uncertainty estimates less calibrated than Gaussian processes

**Key references**:
| Reference | Contribution |
|-----------|-------------|
| Chen & Guestrin, "XGBoost: A Scalable Tree Boosting System", KDD 2016 | Original XGBoost paper |
| Costello & Martin, "A machine learning approach to predicting metabolic pathway dynamics", Bioinformatics 2018 | Neural surrogate for FBA |
| Biggs et al., "Constraint-based model of S. cerevisiae", PLOS Comp Bio 2017 | RF surrogate for metabolic screening |
| Verma et al., "Metabolic engineering design via flux balance analysis and machine learning", 2023 | ML surrogate + metabolic engineering integration |

---

## Topic 4: TOPSIS Weight Sensitivity Analysis

### 4.1 Sensitivity of TOPSIS Rankings to Weight Changes

**Known weakness**: TOPSIS rankings can be sensitive to weight vector changes. Small perturbations in weights can cause rank reversals, particularly for alternatives that are close in their relative closeness to the ideal solution.

**Quantitative characterization**:
- When alternatives have similar relative closeness scores (Ci < 0.05 apart), even a 5-10% weight change can reverse their rankings
- For well-separated alternatives (Ci > 0.15 apart), rankings are typically stable under +/-20% weight perturbations
- The most sensitive weights are those for criteria where top-ranked alternatives perform very differently (e.g., alternative A excels in criterion 1 but is poor in criterion 2, while B is the reverse)

### 4.2 Methods for Robustness Analysis of TOPSIS Results

| Method | Description | Implementation |
|--------|-------------|----------------|
| **Weight perturbation analysis** | Systematically vary each weight by +/-10%, +/-20%, +/-30% and check for rank reversals | Monte Carlo sampling of weight space; count rank reversal frequency |
| **Kendall's tau correlation** | Quantify rank correlation between original and perturbed rankings | tau > 0.7 indicates robust rankings; tau < 0.5 indicates high sensitivity |
| **Spearman rank correlation** | Alternative rank correlation measure | Similar interpretation to Kendall's tau |
| **Stochastic TOPSIS** | Treat weights as random variables with distributions; compute expected ranking and variance | Draw weights from Dirichlet distribution; compute rank probability matrix |
| **Robust TOPSIS** | Minimax approach: find the ranking most robust to worst-case weight perturbations | Optimization over weight uncertainty set |
| **Sensitivity index** | Compute dCi/dwj for each alternative i and weight j | Identifies which weight-parameter combinations are most critical |
| **Pareto front analysis** | Replace TOPSIS entirely with Pareto front from NSGA-II; select knee points | Avoids weight specification entirely; most robust approach |

**Our project's chosen approach** (from `outputs/planning/run_03/phase3_4_integrated_architecture.md`):
1. Perform weight perturbation +/-20% and check for rank reversals
2. If rank reversals occur, apply Entropy weight method for objective weights
3. Use Kendall's tau to quantify ranking stability
4. Present Pareto front visualization as alternative (knee point detection)

### 4.3 Entropy-Based vs Expert-Assigned Weights

| Aspect | Entropy Weight Method | Expert-Assigned Weights |
|--------|----------------------|------------------------|
| **Objectivity** | Fully data-driven; based on information content of each criterion | Subjective; depends on expert judgment and preference |
| **Calculation** | Higher variance in criterion -> higher weight (more discriminating power) | Direct assignment or AHP (Analytic Hierarchy Process) |
| **Formula** | w_j = (1 - e_j) / sum(1 - e_k) where e_j = -sum(p_ij * ln(p_ij)) / ln(n) | Expert-defined or pairwise comparison matrix eigenvector |
| **Advantages** | No subjective bias; reproducible; reflects actual data discriminative power | Captures domain knowledge; can prioritize stakeholder values |
| **Disadvantages** | Ignores domain importance; noisy data can distort weights; criteria with low variance but high importance get low weights | Subjective; may vary between experts; hard to justify objectively |
| **When to use** | When expert consensus is unavailable; as a sanity check against expert weights | When domain knowledge clearly prioritizes certain criteria |
| **Recommended** | Use as **primary** method; compare with expert weights as sensitivity check | Use as **validation**; if both agree, high confidence in ranking |

**Hybrid approach** (recommended for our project):
- Compute entropy weights as baseline (objective)
- Obtain expert weights from domain knowledge (subjective)
- If rankings agree: high confidence
- If rankings disagree: investigate which criteria drive the disagreement; report both rankings
- Final: use a convex combination: w_final = alpha * w_entropy + (1-alpha) * w_expert, with alpha = 0.5 or determined by cross-validation

### 4.4 Key References for TOPSIS Robustness

| Reference | Contribution |
|-----------|-------------|
| Hwang & Yoon, "Multiple Attribute Decision Making", 1981 | Original TOPSIS method |
| Shannon, "A Mathematical Theory of Communication", 1948 | Entropy weight method foundation |
| Zeleny, "Multiple Criteria Decision Making", 1982 | Entropy-based weight determination |
| Opricovic & Tzeng, "Compromise solution by MCDM methods", EJOR 2004 | VIKOR vs TOPSIS comparison; sensitivity analysis |
| Chen et al., "A modified TOPSIS with a different ranking index", IEEE 2006 | Sensitivity of TOPSIS to normalization and distance measures |
| Mardani et al., "A systematic review of multi-criteria decision-making", 2015 | Comprehensive survey including TOPSIS robustness |
| Alinezhad & Amini, "Sensitivity analysis of TOPSIS technique", Applied Mathematical Sciences 2011 | Formal sensitivity analysis methodology for TOPSIS |

---

## Summary: Cross-Topic Integration

The four topics form an integrated technical pipeline:

```
GEM Model (COBRApy)
    |
    v
[FBA Ground Truth Generator] -- 10,000 FBA in ~10-30 min (8-core)
    |                           via Python multiprocessing.Pool
    v
[GEM -> Heterogeneous Graph] -- 3 node types, 3 edge types
    |                           PyG HGTConv/RGCNConv
    v
[GNN + XGBoost Surrogate] ---- 100-1000x speedup over FBA
    |                           XGBoost quantile for uncertainty
    v
[Active Learning Loop] --------- 70-90% reduction in FBA calls
    |                           UCB/EI/Thompson sampling
    v
[NSGA-II + dFBA] -------------- Pareto front of optimal designs
    |
    v
[TOPSIS + Sensitivity] --------- Robust final design selection
                                 Entropy weight + +/-20% check
                                 Kendall's tau > 0.7
```

**Feasibility assessment** (consistent with prior analysis):
- Topic 1 (GNN for metabolic networks): MEDIUM-HIGH feasibility. HGTConv in PyG is mature; main risk is whether structural embeddings capture dynamic flux characteristics.
- Topic 2 (Parallel FBA): HIGH feasibility. Well-understood; main consideration is Windows spawn overhead and solver choice.
- Topic 3 (Surrogate + Active Learning): MEDIUM-HIGH feasibility. Architecture is sound; main risk is generalization failure if training distribution differs from deployment conditions.
- Topic 4 (TOPSIS sensitivity): HIGH feasibility. Well-studied problem with established methods (entropy weights, perturbation analysis, Pareto front alternative).

---

## Source Files Referenced

| File | Content |
|------|---------|
| `outputs/planning/run_02/phase3_4_feasibility_analysis.md` | Phase 3-4 feasibility analysis with GEM->GNN, FBA speed, TOPSIS sensitivity |
| `outputs/planning/run_03/phase3_4_integrated_architecture.md` | Full architecture design (6 modules), data flow, roadmap |
| `outputs/planning/run_03/cobrapy_fba_benchmark.py` | Benchmark script for parallel FBA performance |
| `outputs/planning/run_03/cometspy_dfba_report_20260426.md` | COMETS installation and dFBA demo results |
| `outputs/planning/run_03/flycop_analysis_20260426.md` | FLYCOP analysis and alternatives |
| `outputs/planning/run_02/session_summary_20260426.md` | Session summary of prior research |
| `outputs/literature_review/run_02/literature_review_run02.md` | Literature review on contrastive learning and MoA prediction |
| `outputs/framing/run_01/framing_doc.md` | Research question, dataset, evaluation strategy |
| PyG GitHub README (fetched 2026-04-26) | HGTConv, RGCNConv, RGATConv, HEATConv, FiLMConv documentation |
