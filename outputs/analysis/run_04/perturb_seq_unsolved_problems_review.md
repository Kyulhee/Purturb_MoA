# Core Unsolved Problems in Perturb-seq / Single-Cell Perturbation Biology (2024-2026)

**Date:** 2026-04-27
**Scope:** Literature synthesis of fundamental limitations, computational challenges, and gaps between current methods and biological needs

---

## Executive Summary

Despite rapid methodological advances (20+ new methods in 2025-2026 alone), single-cell perturbation prediction faces 5 deep, interconnected unsolved problems. New methods consistently improve benchmark scores while failing to address structural limitations in evaluation, generalization, and biological fidelity. The field exhibits a pattern of "benchmark chasing" that masks fundamental gaps.

---

## The 5 Core Unsolved Problems

### Problem 1: Evaluation Metrics Are Broken -- We Cannot Distinguish Good Predictions from Artifacts

**The Problem:** Current evaluation metrics for perturbation prediction are incomplete and can be gamed. Methods that score well on benchmarks may actually distort biological structure.

**Evidence:**
- Wang, Leskovec & Regev (Nature Biotechnology, 2026; PMID: 40500472) demonstrated that a trivial 3-layer perceptron ("Islander") **outperforms all leading embedding methods** on standard metrics, yet **distorts biological structures**, making it useless for biological discovery. They introduced scGraph as a corrective metric.
- Csendes et al. (BMC Genomics, 2025; PMID: 40269681) showed that **current Perturb-seq benchmark datasets exhibit low perturbation-specific variance**, making them suboptimal for evaluating models. Many perturbations are indistinguishable from control.
- The Perturbation Discrimination Score (PDS), used in the ARC Virtual Cell Challenge, is **highly sensitive to the choice of distance metric** (L1/L2 vs. cosine) and to the scale of predicted effects (arXiv:2511.16954). Different metrics can rank methods in opposite order.
- The Shesha metric (arXiv:2604.16642) reveals that **effect magnitude alone is insufficient** for evaluating predictions. Two perturbations with identical magnitude can produce qualitatively different outcomes -- one driving cells coherently, another scattering them. Directional coherence (stability) provides a complementary axis that current metrics ignore.
- The GGE framework (arXiv:2603.11233) was proposed specifically to standardize evaluation, indicating the field recognizes this as a crisis.

**What Has Been Tried:**
| Approach | Limitation |
|----------|-----------|
| MSE / R2 on gene expression | Ignores distributional shape; sensitive to highly variable genes |
| Pearson/Spearman correlation | Measures linear association, not distributional fidelity |
| DE overlap / DEG recovery | Binary thresholding loses continuous information |
| PDS (Virtual Cell Challenge) | Sensitive to distance metric choice and effect scaling |
| PDCorr | Improvement over PDS but still incomplete |

**Why It Matters:** Without trustworthy evaluation, the field cannot determine whether any method is actually making progress. New SOTA claims may reflect metric artifacts rather than genuine biological insight.

---

### Problem 2: Combinatorial Explosion with No Sufficient Compositional Theory

**The Problem:** The space of combinatorial perturbations grows exponentially, and no existing method can reliably predict combinatorial outcomes from single-perturbation data alone.

**Scale of the Problem:**
- For K genes, double-knockout space is K(K-1)/2. For iJO1366 (1,367 genes), this is ~934K combinations; triple knockouts exceed 400M.
- Tang et al. (bioRxiv, 2026; PMID: 41648157) note that "probing regulatory elements, interpreting genetic variants, and mapping genetic interactions all challenge the sensitivity and scalability of existing approaches."
- Norman et al. (Science, 2019) showed that **distinct genetic interactions in markedly different parts of the GI manifold may result in similar outcomes when viewed only at the level of fitness** -- high-dimensional phenotyping resolves what simpler readouts cannot.

**Key Non-Additivity Findings:**
- CBL/CNN1 synergy: Individually, neither gene's overexpression strongly induced erythroid markers, but combined perturbation produced a dramatically different transcriptional state (Norman et al., 2019).
- Pleiotropic regulators like CEBPA and GATA1 incur a **"geometric tax"** -- they produce large but incoherent transcriptional shifts, while lineage-specific factors like KLF1 generate tightly coordinated responses (Shesha, arXiv:2604.16642).
- Non-additivity is the rule, not the exception: gene interactions cannot be decomposed into independent contributions.

**What Has Been Tried:**
| Method | Approach | Key Limitation |
|--------|----------|---------------|
| GEARS (Roohani et al., Nature Biotech, 2023) | Graph neural network with gene regulatory network | **Cannot reliably predict combinatorial outcomes from single-gene data alone**; some combinatorial training data is essential; not designed for cross-cell-type transfer |
| CPA (Lotfollahi et al., Molecular Systems Biology, 2023) | Compositional perturbation autoencoder | When Panobinostat single-drug observations were held out, **median DEG score drops from 0.85 to 0.38**; performance depends heavily on training data coverage |
| SAMS-VAE (arXiv:2311.02794) | Sparse additive mechanism shift | Additive decomposition assumption fails for synergistic interactions; sparse latent variables may not capture non-linear epistasis |
| scDFM (arXiv:2602.07103) | Differential flow matching with gene interaction graphs | Reduces MSE by 19.6% in combinatorial settings but still requires some combinatorial training data |
| NAIAD (arXiv:2411.12010) | Active learning for combinatorial discovery | Outperforms by 40% but is an **active learning loop**, not a zero-shot combinatorial predictor; still requires experimental iterations |
| scBIG (arXiv:2602.04901) | Module-inductive representations with flow matching | 6.7% average improvement, but gene program clustering is unsupervised and may not reflect true biological modularity |

**Why It Matters:** Biology needs predictions for perturbations that have never been measured. The combinatorial space is so vast that exhaustive experimentation is impossible. Current methods still require at least some combinatorial training data, making them interpolation tools rather than extrapolation tools.

---

### Problem 3: Out-of-Distribution Generalization and Cross-Cell-Type Transfer

**The Problem:** Models trained on perturbations in one cell type or context fail when applied to different cell types, species, or experimental conditions. This is the single largest barrier to practical utility.

**Evidence:**
- GEARS explicitly states it is **"not designed to handle training across multiple cell types or cross-cell type transfer of predictions"** (GitHub README).
- Sivakumar et al. (Stem Cell Reports, 2025; PMID: 41237780) found that **technical challenges have limited Perturb-seq application in stem-cell-based systems**, highlighting the cell-type dependency problem.
- The Cell-JEPA paper (arXiv:2602.02093) showed that foundation model pretraining **improves absolute-state reconstruction but not effect-size estimation**, suggesting that transfer helps represent cells but not perturbation effects.
- Intermediate-layer analysis (arXiv:2604.14838) found that **optimal embedding layers are task- and context-dependent** -- perturbation optima shift "0-96% across T cell activation states," and first-layer embeddings outperform all deeper layers in quiescent cells, challenging the assumption that deeper layers capture more biologically meaningful features.
- The SynthPert paper (arXiv:2509.25346) achieved **87% accuracy on unseen RPE1 cells** using synthetic reasoning traces, but this was for a simpler question-answering task, not full expression prediction.

**What Has Been Tried:**
| Method | Approach | Key Limitation |
|--------|----------|---------------|
| CellOT (Mittal et al., Nature Methods, 2023) | Neural optimal transport between control/perturbed distributions | Cross-patient generalization demonstrated but cross-cell-type transfer is limited; unpaired distributions require re-training |
| trVAE (arXiv:1910.01791) | Conditional VAE with MMD regularization | Improved Pearson correlations (0.89->0.97) but only for within-dataset conditions |
| TxPert (arXiv:2505.14919) | Knowledge graphs for OOD prediction | Uses biological knowledge networks to generalize, but knowledge graph completeness limits coverage |
| SP-FM (arXiv:2601.11827) | Condition-dependent base distributions for flow matching | Improves OOD for unseen perturbations but not explicitly for unseen cell types |
| C3TL (arXiv:2503.10171) | Causal inductive biases + bulk molecular data | Competitive with large foundation models but still relies on shared causal structure across contexts |
| STRAND (arXiv:2502.06027) | Regulatory DNA sequence conditioning | Expands genomic coverage from ~1.5% to ~95% by conditioning on DNA sequence rather than gene identifiers -- a promising direction for cross-context generalization |
| Mix-Geneformer (arXiv:2507.07454) | Cross-species (human/mouse) unified model | Matches SOTA but variability in zero-shot transfer remains |

**Why It Matters:** In drug discovery and therapeutic development, we need to predict how perturbations will behave in the specific patient cell types we care about -- which are never the cell types used in screens. Without cross-cell-type transfer, Perturb-seq remains a tool for basic science rather than translational applications.

---

### Problem 4: Perturbation Prediction Is Fundamentally Hard -- and We Lack Theory for Why

**The Problem:** The reasons perturbation prediction is hard go beyond data scarcity. The underlying biology exhibits properties that make prediction inherently difficult, and we lack theoretical frameworks to characterize when prediction is even possible.

**Sources of Fundamental Difficulty:**

1. **Non-additivity / Epistasis**: Perturbation effects do not combine linearly. Norman et al. (2019) showed that double-knockout phenotypes cannot be decomposed into independent single-knockout effects. The "GI manifold" reveals diverse interaction types (suppression, synergy, buffering) that require fundamentally different predictive models.

2. **Context-dependence**: The same perturbation produces different effects depending on cell state, cell type, and environmental conditions. The Shesha metric (arXiv:2604.16642) showed that even within the same perturbation, individual cells can move in incoherent directions -- meaning the "effect of a perturbation" is not a single vector but a distribution over possible responses.

3. **High-dimensional output with sparse signal**: Each cell's response is a vector over ~18,000-20,000 genes, but most perturbations affect only a small fraction. The signal-to-noise ratio is extremely low. Nadig et al. (Nat Genet, 2025; PMID: 40259084) explicitly notes that Perturb-seq data "are noisy, and many effects may go undetected."

4. **Pleiotropy and geometric incoherence**: Master regulators like CEBPA and GATA1 produce large but directionally incoherent shifts (the "geometric tax" from Shesha). This means even successful predictions of effect magnitude may miss the critical directional structure.

5. **Distributional heterogeneity**: Perturbations do not shift cell populations as rigid bodies. PerturbDiff (arXiv:2602.19685) showed that modeling distributions in Hilbert space captures "population-level response shifts across hidden factors" that cell-level models miss. The response is not just "where cells go" but "how the shape of the distribution changes."

**What We Do Not Have:**
- No information-theoretic bounds on perturbation predictability
- No theory connecting gene regulatory network structure to perturbation predictability
- No characterization of which perturbation effects are inherently unpredictable vs. merely hard to measure
- No framework for when additive approximations are sufficient vs. when higher-order interactions dominate

---

### Problem 5: The Gap Between What Perturb-seq Measures and What Biology Needs

**The Problem:** There is a fundamental mismatch between the perturbation types, readouts, and scales available in Perturb-seq experiments and the predictions biology actually requires.

**Measurement Gaps:**

| What Perturb-seq Measures | What Biology Needs | Gap |
|---------------------------|-------------------|-----|
| Gene knockouts (CRISPRi/CRISPRko) in immortalized cell lines | Drug responses in primary patient cells | Different perturbation modalities, different cellular contexts |
| Single time-point snapshots (endpoint) | Temporal trajectories of response | Perturb-seq captures state, not dynamics; dFBA/simulation needed for trajectories |
| Transcriptome only (scRNA-seq) | Multi-omic response (proteome, epigenome, metabolome) | Single-modality measurement misses most of the biological response |
| 1-2 gene perturbations | Multi-gene programs, pathway-level perturbations | Combinatorial coverage is minuscule relative to the space |
| Bulk-level effects of perturbations | Cell-type-specific and cell-state-specific effects | Population averaging obscures heterogeneity |
| In vitro cell lines | In vivo tissue context with microenvironment | Tissue context fundamentally changes perturbation responses |

**Evidence:**
- Zhang et al. (BMC Genomics, 2026; PMID: 41826830) found that pooled CRISPRi screens have "broader application limited by technical" factors.
- Replogle et al. (Cell, 2022; PMID: 35688146) scaled Perturb-seq to genome-wide levels but acknowledged that "high-content phenotypic screens... have been used at limited scales."
- Anglada-Girotto et al. (Nucleic Acids Res, 2025; PMID: 41101775) note that "systematically mapping disease-driver regulatory interactions at large scale remains challenging."
- Tang et al. (2026) highlight that "probing regulatory elements, interpreting genetic variants, and mapping genetic interactions" all push beyond what current screen sensitivity can deliver.

**The Drug Response Problem:**
- scPPDM (arXiv:2510.11726) is "the first diffusion-based framework for single-cell drug-response prediction," published in late 2025, indicating this gap persisted for years.
- CPA can incorporate chemical representations to predict "response to completely unseen drugs," but performance drops from 0.85 to 0.38 when a drug is held out.
- The Tahoe-100M dataset (used by SCALE, scPPDM, and others) represents the largest drug perturbation benchmark, but covers only a fraction of chemical space.

---

## Comprehensive Method Landscape (2023-2026)

### Gene Perturbation Prediction Methods

| Method | Year | Architecture | Key Innovation | Persistent Limitation |
|--------|------|-------------|----------------|----------------------|
| GEARS | 2023 | GNN + gene regulatory network | Graph-enhanced combinatorial prediction | No cross-cell-type transfer; needs combinatorial training data |
| CPA | 2023 | Compositional autoencoder | Modular architecture for unseen combinations | Severe degradation for unseen drugs (0.85->0.38) |
| CellOT | 2023 | Neural optimal transport | Unpaired distribution mapping | Requires re-training per perturbation; limited cross-context transfer |
| SAMS-VAE | 2023 | Sparse additive VAE | Disentangled perturbation-specific latents | Additive assumption fails for strong epistasis |
| PerturbNet | 2023 | Flow-based + GNN | High R2 on sci-Plex | Limited generalization evidence |
| scDFM | 2026 | Differential flow matching | Perturbation-Aware Differential Transformer | 19.6% MSE reduction but still needs combinatorial training data |
| scBIG | 2026 | Module-inductive + flow matching | Gene program clustering | Unsupervised clustering may not reflect biological modularity |
| PerturbDiff | 2026 | Distributional diffusion in Hilbert space | Models full population-level distributions | Computationally expensive; limited to datasets with many cells per perturbation |
| PRiMeFlow | 2026 | End-to-end flow matching | Won ARC Virtual Cell Challenge Generalist Prize | Challenge benchmarks have known metric issues (PDS sensitivity) |
| Lingshu-Cell | 2026 | Masked discrete diffusion | Discrete token space for sparse data | Leading on VCC H1 but evaluation metric concerns persist |
| SCALE | 2026 | Conditional transport + LLaMA encoding | 12.51x pretrain speedup; SOTA on Tahoe-100M | Still evaluated on metrics with known flaws |
| Departures | 2025 | Neural Schrodinger bridges | Minibatch-OT pairing for distribution alignment | SOTA but Schrodinger bridge optimization is unstable |
| PRESCRIBE | 2025 | Deep evidential regression | Epistemic + aleatoric uncertainty quantification | Uncertainty is calibrated but does not fix prediction errors |

### Drug / Chemical Perturbation Methods

| Method | Year | Architecture | Key Innovation | Persistent Limitation |
|--------|------|-------------|----------------|----------------------|
| scPPDM | 2025 | Diffusion + factorized guidance | First diffusion for drug response; SOTA on Tahoe-100M | Only evaluated on in vitro datasets |
| CellFlux | 2025 | Flow matching for morphology | Distribution-wise transformations | Morphology prediction, not transcriptomic |
| MolPhenix | 2024 | Contrastive learning (molecule + phenotype) | 8.1x improvement in zero-shot molecular retrieval | Retrieval, not de novo prediction |

### Cross-Context / Foundation Model Approaches

| Method | Year | Key Innovation | Persistent Limitation |
|--------|------|----------------|----------------------|
| scGPT (perturbation) | 2024 | Foundation model zero-shot | Attention provides "no incremental value for perturbation prediction" over trivial baselines (arXiv:2602.17532) |
| Cell-JEPA | 2026 | Joint-embedding predictive architecture | 36% improvement in zero-shot transfer (AvgBIO) but does NOT improve effect-size estimation |
| STRAND | 2026 | DNA sequence conditioning | Expands coverage from 1.5% to 95% but sequence->expression mapping is still imperfect |
| C3TL | 2026 | Causal inductive biases + bulk data | Competitive with foundation models but assumes shared causal structure |
| TxPert | 2025 | Knowledge graph-enhanced OOD | Depends on knowledge graph completeness |
| Mix-Geneformer | 2025 | Cross-species unified model | Variability in zero-shot transfer |
| GenoHoption | 2024 | Gene network graphs + foundation model | 3.86% improvement on perturbation -- marginal |

### LLM / Agent Approaches

| Method | Year | Key Innovation | Persistent Limitation |
|--------|------|----------------|----------------------|
| PBio-Agent | 2026 | Multi-agent LLM with difficulty-aware sequencing | Improves PerturbQA but QA != expression prediction |
| SynthPert | 2025 | Synthetic reasoning traces for LLM fine-tuning | 87% on unseen cell types for QA, not expression prediction |
| VCWorld | 2025 | White-box simulator + LLM reasoning | Interpretable but data-efficient only for simple cascades |
| CellForge | 2025 | Multi-agent architecture design | Automated but limited to existing architectural primitives |

---

## Key Insight: The Evaluation Crisis Undermines All Progress Claims

The most meta-level finding from this review is that the field has a **measurement problem before it has a prediction problem**:

1. **Low perturbation-specific variance** in benchmark datasets means many perturbations look like control (Csendes et al., 2025). Models trained on such data cannot be meaningfully evaluated.
2. **Metric sensitivity**: PDS ranks methods differently depending on distance metric choice (arXiv:2511.16954). A method can be SOTA under L2 but below-baseline under cosine.
3. **Good metrics /= biological validity**: Islander outperforms all methods on standard metrics while distorting biological structure (Wang et al., Nature Biotech, 2026).
4. **Magnitude vs. coherence**: The Shesha metric shows that effect magnitude and directional coherence are distinct axes. Models optimizing only for magnitude accuracy may produce biologically meaningless predictions.

This means the apparent rapid progress (20+ new methods in 2025-2026) may be partially illusory. Without fixing evaluation, the field risks a "AI benchmarking trap" where methods improve on flawed metrics without advancing biological understanding.

---

## Papers Referenced (25+ papers)

1. Wang H, Leskovec J, Regev A. "Limitations of cell embedding metrics assessed using drifting islands." Nature Biotechnology, 2026. PMID: 40500472
2. Csendes G, Sanz G, Szalay KZ, Szalai B. "Benchmarking foundation cell models for post-perturbation RNA-seq prediction." BMC Genomics, 2025. PMID: 40269681
3. Tang A, Ardy RC, Mendes RE, Norman TM. "Scaling perturbations: beyond genome-scale CRISPR screens." bioRxiv, 2026. PMID: 41648157
4. Zhang H et al. "Insights from pooled CRISPRi single-cell screens in K562 cells." BMC Genomics, 2026. PMID: 41826830
5. Sivakumar S et al. "Benchmarking and optimizing Perturb-seq in differentiating human pluripotent stem cells." Stem Cell Reports, 2025. PMID: 41237780
6. Roohani Y, Huang K, Leskovec J. "Predicting transcriptional outcomes of novel multigene perturbations with GEARS." Nature Biotechnology, 2023.
7. Lotfollahi M et al. "Predicting cellular responses to complex perturbations in high-throughput screens." Molecular Systems Biology, 2023.
8. Mittal S et al. "Learning single-cell perturbation responses using neural optimal transport." Nature Methods, 2023.
9. Norman TM et al. "Combinatorial CRISPR screens reveal genetic interactions in human cells." Science, 2019.
10. Replogle JM et al. "Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq." Cell, 2022. PMID: 35688146
11. Nadig R et al. "TRADE: transcriptome-wide impact estimation." Nature Genetics, 2025. PMID: 40259084
12. Anglada-Girotto M et al. "Using single-cell perturbation screens to decode the regulatory architecture of splicing factor programs." Nucleic Acids Research, 2025. PMID: 41101775
13. SAMS-VAE. arXiv:2311.02794 (2023)
14. scDFM. arXiv:2602.07103 (2026) - ICLR 2026
15. scBIG. arXiv:2602.04901 (2026)
16. NAIAD. arXiv:2411.12010 (2024)
17. PerturbDiff. arXiv:2602.19685 (2026)
18. PRiMeFlow. arXiv:2604.13986 (2026)
19. Lingshu-Cell. arXiv:2603.25240 (2026)
20. SCALE. arXiv:2603.17380 (2026)
21. Departures. arXiv:2511.11265 (2025)
22. PRESCRIBE. NeurIPS 2025
23. Cell-JEPA. arXiv:2602.02093 (2026)
24. STRAND. arXiv:2502.06027 (2026)
25. C3TL. arXiv:2503.10171 (2026)
26. TxPert. arXiv:2505.14919 (2025)
27. SP-FM. arXiv:2601.11827 (2026)
28. Shesha (geometric coherence). arXiv:2604.16642 (2026)
29. PDS sensitivity analysis. arXiv:2511.16954 (2025)
30. scGPT interpretability evaluation. arXiv:2602.17532 (2026)
31. SynthPert. arXiv:2509.25346 (2025)
32. PBio-Agent. arXiv:2602.07408 (2026)
33. VCWorld. arXiv:2512.00306 (2025) - ICLR 2026
34. scPPDM. arXiv:2510.11726 (2025)
35. SAVE. arXiv:2604.16776 (2026) - ICLR 2026
36. Intermediate layers. arXiv:2604.14838 (2026) - ICLR 2026 Workshop
37. PerturbQA benchmark. arXiv:2502.21290 (2025)
38. LLM survey for virtual cell. arXiv:2510.07706 (2025)
39. How to build virtual cell with AI. arXiv:2409.11654 (2024)
40. GGE evaluation framework. arXiv:2603.11233 (2026)
41. Mix-Geneformer. arXiv:2507.07454 (2025)
42. GenoHoption. arXiv:2411.06331 (2024)
43. trVAE. arXiv:1910.01791 (2019)

---

## Summary: Where the Gap Is Largest

The biggest gap between current methods and what biology needs is **not** in architecture design (new transformers, flow matching, diffusion models) but in three structural problems:

1. **Evaluation validity**: We cannot trust our progress measurements. The metrics are gamed, the benchmarks have low signal, and good metric performance does not imply biological fidelity.

2. **Compositional theory**: We have no theory predicting when perturbation effects compose additively vs. when they exhibit epistasis. Without this, combinatorial prediction remains empirical interpolation, not principled extrapolation.

3. **Context transfer**: We cannot predict how a perturbation learned in one cell type will behave in another. This is the single biggest barrier to translational impact. The STRAND approach (DNA sequence conditioning) and C3TL (causal transfer) are the most promising directions but remain early-stage.

These three problems are deeply interconnected: better evaluation would reveal which methods actually transfer; better compositional theory would constrain the search space for transfer learning; and better transfer would enable training on abundant data (cell lines) to predict in scarce contexts (patient cells).
