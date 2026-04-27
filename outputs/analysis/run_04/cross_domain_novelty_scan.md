# Cross-Domain Novelty Scan: High-Novelty Hypothesis Candidates for Perturb-seq

**Date**: 2026-04-27
**Scope**: 10 cross-domain intersections with single-cell perturbation prediction
**Method**: Systematic arXiv search (2024-2026), key paper deep-reads, gap analysis

---

## Executive Summary

Ten cross-domain ideas were scanned. Three emerged as **high-novelty, under-explored** candidates where no direct Perturb-seq implementation exists but theoretical foundations are mature enough to attempt. Two are **medium-novelty** (early work exists but major gaps remain). Five are **partially addressed** by recent work but with clear next-step opportunities.

| # | Domain | Novelty | Feasibility | Impact | Status |
|---|--------|---------|-------------|--------|--------|
| 1 | Causal Representation Learning | HIGH | Medium | Very High | FCR exists but disentangled *causal* perturbation representations are untested |
| 2 | In-Context Learning / Foundation Models | MEDIUM | High | Very High | MapPFN (Jan 2026) exists; open problems remain |
| 3 | Meta-Learning / Few-Shot | HIGH | Medium | High | No direct Perturb-seq meta-learning work found |
| 4 | Equivariant Neural Networks | HIGH | Low-Medium | Medium | Zero work at intersection of equivariance + perturbation prediction |
| 5 | Neural ODE / CNF | MEDIUM | Medium | High | Cell-MNN (Oct 2025) applies ODE latent; perturbation-specific ODE untested |
| 6 | Mechanistic Interpretability | HIGH | Medium | High | Zero work interpreting what perturbation models learn |
| 7 | Optimal Transport | LOW | High | Medium | Saturated; W1 solver, Conditional Monge Gap, SCALE all exist |
| 8 | Compositional Generalization | HIGH | Medium | Very High | No principled compositional framework for combinatorial perturbations |
| 9 | Physics-Informed / Biologically-Constrained NNs | HIGH | Medium | High | No work embedding GRN/pathway topology as hard constraints |
| 10 | Invariant Risk Minimization / OOD | MEDIUM | Medium | High | Cross-cell-type work emerging (CFM-GP, C3TL) but no IRM/Causal approach |

---

## 1. Causal Representation Learning

### What exists
- **FCR** (Mao et al., Oct 2024, arXiv:2410.22472): Factorized Causal Representations that disentangle z_x (covariate), z_t (treatment), z_tx (interaction) with identifiability proofs. Applied to single-cell perturbation data across cell lines. **Closest to the idea.**
- **Uhler & Zhang** (Nov 2025, arXiv:2511.04790): ICM 2026 survey outlining framework for causal structure + representation learning from multi-modal perturbation data. Defines three problems: causal discovery, learning causal variables, optimal perturbation design.
- **Reizinger et al. (Schölkopf group)** (Jun 2024, arXiv:2406.14302): IEM framework unifying causal structure and representation learning via exchangeable mechanisms. Relaxed identifiability conditions. ICLR 2025.
- **Rajendran et al. (Schölkopf group)** (Feb 2024, arXiv:2402.09236): Unifies causal representation learning with foundation models. Provably recoverable "concepts" from diverse data. NeurIPS 2024.

### What is novel / untested
1. **Causal interventions on the latent space itself**: FCR learns disentangled factors but does not perform *do-calculus* on them. What if we could intervene on z_tx directly to predict combinatorial effects?
2. **Independent Causal Mechanisms (ICM) as inductive bias for perturbation transfer**: The Schölkopf group's ICM principle (mechanisms are autonomous/invariant) has never been used to regularize perturbation prediction across cell types. If the mechanism z_tx is truly causal, it should be invariant across cell contexts -- enabling zero-shot transfer.
3. **Causal discovery from perturbation embeddings**: Can we reverse-engineer the causal graph (which genes regulate which) from the learned treatment representations z_t? This is the "learning causal variables" problem from Uhler & Zhang, applied specifically to perturbation prediction models.

### Assessment
- **Novelty**: HIGH. FCR is the first step, but no one has connected ICM invariance to cross-cell-type transfer, or used perturbation model latents for causal discovery.
- **Feasibility**: MEDIUM. Requires identifiability guarantees that are hard to satisfy with real data. FCR's block-wise identifiability is a start but strong assumptions needed.
- **Impact**: VERY HIGH. If z_tx is truly invariant, zero-shot cross-cell-type prediction becomes possible without any target-cell measurements.
- **Core hypothesis**: "The interaction representation z_tx learned by FCR is invariant across cell types by the ICM principle; therefore, z_tx learned from cell type A can predict perturbation effects in cell type B without any fine-tuning."

---

## 2. In-Context Learning / Foundation Models

### What exists
- **MapPFN** (Sextro et al., Jan 2026, arXiv:2601.21092): Prior-data fitted network that uses in-context learning to predict post-perturbation distributions. Pre-trained on synthetic causal priors (in silico gene knockouts). Zero-shot performance on par with models trained on real data. **This IS the idea, implemented.**
- **PT-RAG** (Di Francesco et al., Mar 2026, arXiv:2603.07233): Retrieval-augmented generation for perturbation prediction. Differentiable, cell-type-aware retrieval via Gumbel-Softmax. Key finding: vanilla RAG *fails* -- differentiable retrieval is essential.
- **SynthPert** (Phillips et al., Sep 2025, arXiv:2509.25346): Synthetic reasoning distillation for LLM-based perturbation prediction. 87% accuracy on unseen cell type (RPE1).
- **Civale et al.** (Apr 2026, arXiv:2604.14838): Intermediate layers of scFoundation models outperform final layers. First-layer embeddings beat all deeper layers in quiescent cells. Layer selection matters.

### What is novel / untested
1. **Scaling laws for perturbation PFNs**: MapPFN was demonstrated on small datasets. What are the scaling properties? Does in-context perturbation prediction follow Chinchilla-like scaling laws?
2. **Multi-modal in-context perturbation**: No work combines gene expression + protein + imaging perturbation data in a single in-context framework.
3. **Active learning via in-context evidence**: MapPFN can incorporate new evidence at inference time. This enables closed-loop experimental design -- design experiment, observe result, feed into context, predict next experiment. No one has built this loop.

### Assessment
- **Novelty**: MEDIUM (MapPFN already implements the core idea).
- **Feasibility**: HIGH (transformer architectures well-understood, synthetic data generation feasible).
- **Impact**: VERY HIGH if scaling works.
- **Core hypothesis**: "A PFN pre-trained on a sufficiently rich prior over causal perturbation mechanisms can predict combinatorial perturbation effects zero-shot by composing individual perturbation contexts in-context."

---

## 3. Meta-Learning / Few-Shot Perturbation Prediction

### What exists
- **MetaCaDI** (Ong et al., Oct 2025, arXiv:2510.22298): Meta-learning for causal discovery with few-shot intervention targets. Recovers causal graphs from as few as 10 data instances. Applied to gene expression. **Closest existing work.**
- MapPFN (above) can be seen as a form of meta-learning via in-context learning.

### What is novel / untested
1. **MAML/ProtoNet-style meta-learning for perturbation effects**: No work treats perturbation prediction as a meta-learning problem where each gene knockout is a "task" and the goal is to learn to predict new knockouts from 5-10 cells. This is a direct application of few-shot learning that is completely unexplored.
2. **Task augmentation via synthetic perturbation data**: Generate diverse perturbation tasks synthetically (e.g., from ODE gene regulatory models) to augment the meta-training distribution, similar to how MapPFN uses synthetic priors but with explicit task-level structure.
3. **Cross-species meta-transfer**: Meta-learn across organisms (human, mouse, zebrafish perturbation data) so that new perturbations in any species require minimal measurements.

### Assessment
- **Novelty**: HIGH. No MAML/ProtoNet-style few-shot learning for perturbation prediction exists.
- **Feasibility**: MEDIUM. Requires sufficient task diversity for meta-training. Perturb-seq datasets cover ~100-1000 genes, which may be enough.
- **Impact**: HIGH. Reducing measurement requirements by 10-100x would transform experimental design.
- **Core hypothesis**: "Perturbation effects share learnable structure across genes (e.g., pathway co-membership, transcription factor co-regulation); a meta-learner can extract this structure and predict new gene perturbation effects from k=5-10 post-perturbation cells."

---

## 4. Equivariant Neural Networks

### What exists
- Equivariant GNNs are mature in molecular modeling (E2Former, TransIP, EGNN-based methods) but have **zero application** to single-cell perturbation prediction.
- No work on identifying biological symmetries in perturbation response that should be respected by neural architectures.

### What is novel / untested
1. **Permutation-equivariance for gene program structure**: Genes within a regulatory program can be reordered without changing the program-level effect. A network that is equivariant to within-program permutations would naturally capture the modular structure that scBIG (arXiv:2602.04901) explicitly models.
2. **Scale-equivariance for dose-response**: Drug perturbations at different doses produce scaled responses. A scale-equivariant architecture would naturally extrapolate dose-response without seeing all doses.
3. **Equivariance to known regulatory ordering**: If gene A activates gene B, perturbation of A should produce effects that are consistent with (and supersede) perturbation of B alone. This is a form of partial-order equivariance -- entirely unexplored.

### Assessment
- **Novelty**: VERY HIGH. Zero intersection between equivariant NNs and perturbation prediction.
- **Feasibility**: LOW-MEDIUM. Identifying the correct symmetry group for biological perturbations is non-trivial. The partial-order equivariance idea is theoretically appealing but hard to implement.
- **Impact**: MEDIUM. Equivariance would improve data efficiency and generalization, but the magnitude of improvement is uncertain.
- **Core hypothesis**: "A neural network that is equivariant to within-gene-program permutations will learn more data-efficient representations of perturbation effects and generalize better to combinatorial perturbations where the same program is perturbed in different ways."

---

## 5. Neural ODE / Continuous Normalizing Flows

### What exists
- **Cell-MNN** (von Bassewitz et al., Oct 2025, arXiv:2510.02903): Encoder-decoder where latent space is a locally linearized ODE governing cellular evolution. Learns interpretable gene interactions. ICLR 2026. **Applied to trajectory inference, not perturbation prediction directly.**
- **scDFM** (Yu et al., Feb 2026, arXiv:2602.07103): Distributional flow matching for perturbation prediction. ICLR 2026. Models full distribution, not individual cell dynamics.
- **CFM-GP** (Abir et al., Aug 2025, arXiv:2508.08312): Conditional flow matching for cross-cell-type perturbation prediction.
- **SP-FM** (Rubbi et al., Jan 2026, arXiv:2601.11827): Shortest-path flow matching with mixture-conditioned bases for OOD generalization.
- **CaLMFlow** (He et al., Oct 2024, arXiv:2410.05292): Reformulates flow matching as Volterra integral equation solved by causal language models.

### What is novel / untested
1. **Neural ODE with perturbation-specific vector fields**: Cell-MNN models differentiation dynamics. What if perturbation alters the vector field itself? A Neural ODE where the perturbation modulates the ODE right-hand-side (f(x, t, perturbation)) would model perturbation effects as continuous dynamics -- you could "interpolate" between perturbation strengths, predict time-resolved effects, and compose perturbations by composing vector field modifications.
2. **Continuous-time perturbation trajectories**: All current methods predict the endpoint (post-perturbation state). No method predicts the trajectory of how cells transition from pre- to post-perturbation state over time. This requires time-series perturbation data (which exists for some systems).
3. **Composing vector field perturbations**: If perturbation A modifies f to f_A and perturbation B modifies f to f_B, can we predict f_{A+B} = f + delta_A + delta_B? This is a compositional generalization question expressed in ODE language.

### Assessment
- **Novelty**: MEDIUM. Flow matching is well-explored, but the specific formulation of perturbation-as-vector-field-modulation is untested.
- **Feasibility**: MEDIUM. Requires time-series data for training. Norman et al. (2019) Perturb-seq has no time dimension, but newer datasets (e.g., CROP-seq time course) do.
- **Impact**: HIGH. Continuous-time predictions would be a qualitative advance, not just incremental.
- **Core hypothesis**: "Perturbation effects can be modeled as modifications to the ODE vector field governing cellular dynamics; composing perturbations corresponds to composing vector field modifications, enabling prediction of combinatorial effects."

---

## 6. Mechanistic Interpretability

### What exists
- **Le et al.** (Jul 2024, arXiv:2407.10785): Sparse autoencoders (SAEs) applied to pathology foundation models. Found "monosemantic biological concepts" in individual SAE dimensions. **Closest work.**
- **Conan** (May 2025, arXiv:2505.00555): On mechanistic interpretability of neural networks for causality in bio-statistics. Probes internal representations for causal bio-statistical analysis.
- **PerturBench** (Wu et al., Aug 2024, arXiv:2408.10609): Found mode collapse in widely-used perturbation models and that "simpler architectures are generally competitive." But did not investigate *why* models fail or *what* they learn.

### What is novel / untested
1. **SAE-based dissection of perturbation model latents**: Apply sparse autoencoders to perturbation prediction models (CellOT, scDFM, GEARS, etc.) to discover whether they learn biologically meaningful features (pathway activations, TF activities) or statistical shortcuts. This is completely unexplored.
2. **Causal scrubbing for perturbation models**: Apply the causal scrubbing methodology (from AI safety) to determine which input features causally influence perturbation predictions. Do models use the perturbation identity, or do they rely on confounders (batch, cell cycle state)?
3. **Probing for known biology**: Train linear probes on perturbation model latents to test whether known pathway structure (e.g., MAPK cascade, p53 network) is encoded. If models do not encode real biology, their generalization is on shaky ground.

### Assessment
- **Novelty**: VERY HIGH. Zero work on mechanistic interpretability of perturbation prediction models.
- **Feasibility**: MEDIUM. SAE methodology is well-established from LLM interpretability. Perturbation models are smaller than LLMs, making analysis tractable.
- **Impact**: HIGH. Understanding what models learn would guide architecture design and reveal whether current models are learning biology or shortcuts.
- **Core hypothesis**: "Current perturbation prediction models do not encode known pathway structure in their latent representations; they learn statistical shortcuts that work for in-distribution prediction but fail for combinatorial and cross-cell-type generalization."

---

## 7. Optimal Transport

### What exists (SATURATED domain)
- **CellOT** (Bunne et al., 2022): W2-based neural OT. The original.
- **W1 OT solver** (Chen et al., Nov 2024, arXiv:2411.00614): 25-45x speedup over W2, single maximization instead of min-max. ISMB/ECCB 2025.
- **Conditional Monge Gap** (Driessen et al., Apr 2025, arXiv:2504.08328): Conditional OT maps that generalize to unseen drugs/dosages via cross-task learning.
- **SCALE** (Chen et al., Mar 2026, arXiv:2603.17380): Atlas-level conditional transport with BioNeMo infrastructure. 12x speedup.
- **Distribution-Conditioned Transport** (Fishman et al., Mar 2026, arXiv:2603.04736): Lifts autoencoders to distribution space.
- **STRAND** (Fu et al., Feb 2026, arXiv:2602.10156): Sequence-conditioned transport. Encodes genomic locus for zero-shot inference.
- **Unbalanced Monge Maps** (Eyring et al., 2023/2024, arXiv:2311.15100): Handles unbalanced transport.

### Fundamental limitations and how to overcome them
1. **OT assumes unpaired data but ignores biology**: OT finds the cheapest transport between distributions, but biological constraints (known regulatory relationships, pathway topology) are not incorporated. **Overcome by**: Constrained OT where transport plans must respect known biological graph structure.
2. **Compositional OT is undefined**: OT maps for individual perturbations cannot be composed to predict combinatorial effects. **Overcome by**: Perturbation-conditioned vector fields (Neural ODE approach) where composition is vector field addition.
3. **Cell-level transport is ill-posed**: The "same cell" does not exist pre/post perturbation (destructive measurement). OT resolves this mathematically but the solution is not unique. **Overcome by**: PerturbDiff's distribution-level approach (arXiv:2602.19685) operates on distributions directly.

### Assessment
- **Novelty**: LOW. Highly saturated domain.
- **Feasibility**: HIGH.
- **Impact**: MEDIUM. Incremental improvements at this point.
- **Best remaining opportunity**: Constrained OT incorporating biological graph structure.

---

## 8. Compositional Generalization

### What exists
- **scDFM** (Feb 2026, arXiv:2602.07103): 19.6% MSE reduction in combinatorial setting. Uses distributional flow matching. Does not explicitly model composition.
- **scBIG** (Ruan et al., Feb 2026, arXiv:2602.04901): Module-inductive representation that models coordinated gene programs. 6.7% average improvement over baselines. Implicitly helps composition by structuring representations, but does not have a compositional rule.
- **PerturbDiff** (Yuan et al., Feb 2026, arXiv:2602.19685): Distribution-level diffusion for perturbation. Better generalization to unseen perturbations, but no explicit compositional mechanism.

### What is novel / untested
1. **Neural module networks for perturbation composition**: Inspired by cognitive science work on compositional generalization (Andreas et al., 2016), build neural networks where each perturbation is a "module" and combinatorial effects are computed by composing modules. The key question: what is the right composition operation? Addition in latent space? Multiplication? Something else?
2. **Why composition fails biologically**: Gene perturbations can exhibit epistasis, buffering, synthetic lethality, and redundancy. These are fundamentally non-compositional phenomena. A model that assumes additive composition will fail on epistatic interactions. **No work has systematically characterized which types of gene interactions are compositional and which are not.**
3. **Compositional generalization benchmarks**: No dedicated benchmark exists for testing whether models can predict double-knockout effects from single-knockout data. The Replogle dataset has some double perturbations but is not designed for this evaluation.

### Assessment
- **Novelty**: VERY HIGH. No principled compositional framework exists for perturbation prediction.
- **Feasibility**: MEDIUM. The core difficulty is that biology is often non-compositional (epistasis). But even partial compositional structure would be valuable.
- **Impact**: VERY HIGH. Combinatorial perturbation space grows exponentially; even approximate composition would be transformative.
- **Core hypothesis**: "Gene perturbation effects decompose into pathway-level modules that compose additively for genes in independent pathways and interact multiplicatively for genes in the same pathway; a model that learns this structure can predict combinatorial effects from single-gene perturbation data alone."

---

## 9. Physics-Informed / Biologically-Constrained Neural Networks

### What exists
- **PINN for biological systems**: Several papers use physics-informed losses for biological ODEs (glucose-insulin dynamics, adoptive cell therapy, neural population dynamics). **None apply to perturbation prediction.**
- **scBIG** (arXiv:2602.04901): Implicitly captures gene program structure through data-driven clustering. Does not embed hard biological constraints.
- **STRAND** (arXiv:2602.10156): Encodes genomic sequence as a soft prior. Not a hard constraint.

### What is novel / untested
1. **GRN-constrained perturbation prediction**: Train a perturbation prediction model with a physics-informed loss that penalizes predictions violating known gene regulatory network structure. For example: if gene A is known to activate gene B, the predicted effect of knocking out A should include downregulation of B. Current models have no such constraint.
2. **Pathway topology as a hard constraint**: Embed signaling pathway topology (e.g., KEGG, Reactome) as a directed acyclic graph and require that perturbation effects propagate along this graph. This is a direct application of PINN methodology to biological networks.
3. **Thermodynamic constraints on gene expression**: Gene expression levels are bounded (0 to maximum transcription rate). Protein folding and degradation impose additional constraints. None of these are enforced in current models.

### Assessment
- **Novelty**: VERY HIGH. No work embeds biological network constraints as hard constraints in perturbation prediction models.
- **Feasibility**: MEDIUM. Requires high-quality pathway/GRN databases. KEGG and Reactome exist but are incomplete and sometimes incorrect. Soft constraints may be more practical than hard constraints.
- **Impact**: HIGH. Constrained models would be more interpretable, more data-efficient, and more likely to generalize.
- **Core hypothesis**: "A perturbation prediction model constrained by known gene regulatory network structure will (a) require less training data, (b) produce biologically plausible predictions for out-of-distribution perturbations, and (c) have interpretable failure modes tied to incorrect GRN annotations."

---

## 10. Invariant Risk Minimization / OOD Generalization

### What exists
- **CFM-GP** (Aug 2025, arXiv:2508.08312): Cell type-agnostic flow matching. Single model across cell types.
- **C3TL** (Scholkemper & Mukherjee, Mar 2026, arXiv:2603.13051): Lightweight causal context transfer. Competitive with foundation models using only bulk data.
- **Conditional Monge Gap** (Apr 2025, arXiv:2504.08328): Cross-task OT learning that generalizes to unseen drugs.
- **SP-FM** (Jan 2026, arXiv:2601.11827): Mixture-conditioned bases for OOD flow matching.
- **SynthPert** (Sep 2025, arXiv:2509.25346): 87% accuracy on unseen RPE1 cell type via synthetic reasoning distillation.

### What is novel / untested
1. **IRM/Causal invariance for perturbation transfer**: None of the above methods use invariant risk minimization (Arjovsky et al., 2019) or related causal invariance principles. The core idea: if we treat each cell type as an "environment," features that are invariant across environments should be the ones that generalize. This is a direct application of IRM that is completely untested.
2. **Domain adaptation with biological structure**: Standard domain adaptation assumes arbitrary domain shift. In biology, the shift between cell types is structured (they share pathways, differ in expression levels). No work exploits this structure for domain adaptation.
3. **Generalization bounds for perturbation prediction**: No theoretical work on what makes perturbation prediction generalize across cell types. Information-theoretic bounds could guide model design.

### Assessment
- **Novelty**: MEDIUM. Cross-cell-type transfer is being actively explored, but IRM/causal invariance approaches are absent.
- **Feasibility**: MEDIUM. IRM has known optimization difficulties (gradient starvation, failure in nonlinear settings). But the structured nature of biological domains may help.
- **Impact**: HIGH. Principled OOD generalization would provide theoretical guarantees lacking in current empirical approaches.
- **Core hypothesis**: "Perturbation effect features that are invariant across cell types (in the IRM sense) correspond to the true causal mechanism of perturbation action; training models to extract invariant features will produce zero-shot cross-cell-type predictions."

---

## Top 3 High-Novelty Hypothesis Candidates

### H1: Causal Invariance Hypothesis (Domains 1 + 10)
**The interaction representation z_tx from FCR, when constrained by the Independent Causal Mechanism principle, is invariant across cell types and enables zero-shot perturbation transfer.**

- Testable prediction: z_tx learned from K562 cells will predict perturbation effects in RPE1 cells without fine-tuning, and the prediction accuracy will correlate with the degree of ICM invariance in z_tx.
- Required: FCR implementation + IRM training objective + multi-cell-line perturbation data (Norman, Replogle).
- Why now: FCR provides the disentangled representation; IRM provides the invariance objective; the combination is novel.

### H2: Compositional Module Hypothesis (Domains 3 + 8)
**Gene perturbation effects decompose into pathway-level modules that compose via a learnable interaction function; a meta-learner trained on single-gene perturbations can predict combinatorial effects by composing modules.**

- Testable prediction: A model trained only on single-gene knockouts will predict double-knockout effects within R2 > 0.5 for genes in independent pathways and identify epistatic interactions as residuals from the compositional model.
- Required: Module network architecture + meta-learning training on single-gene data + evaluation on double-knockout data (Norman dataset has combinatorial perturbations).
- Why now: scBIG shows gene program structure helps; MapPFN shows in-context learning works; the composition gap is the remaining frontier.

### H3: Mechanistic Interpretability Audit (Domain 6)
**Current perturbation prediction models learn statistical shortcuts, not biological mechanisms; applying sparse autoencoders will reveal that model latents do not encode known pathway structure.**

- Testable prediction: Linear probes on CellOT/scDFM/GEARS latents will fail to recover known pathway annotations (e.g., MAPK, p53 targets) above chance level; SAE features will be monosemantic for batch/confounder variables but polysemantic or absent for pathway variables.
- Required: SAE implementation + pathway gene sets (MSigDB) + trained perturbation models (open-source CellOT, GEARS).
- Why now: SAE methodology is mature from LLM interpretability; Le et al. (2024) demonstrated it works on biological foundation models; perturbation models are small enough for thorough analysis.

---

## Key Papers Referenced

| arXiv ID | Title | Key Contribution |
|----------|-------|-----------------|
| 2410.22472 | FCR: Factorized Causal Representations | Disentangled z_x, z_t, z_tx with identifiability proofs |
| 2601.21092 | MapPFN: In-Context Perturbation Prediction | PFN with in-context learning for perturbation effects |
| 2603.13051 | C3TL: Causal Cellular Context Transfer | Lightweight causal transfer competitive with FMs |
| 2602.19685 | PerturbDiff: Functional Diffusion | Distribution-level diffusion for perturbation |
| 2602.07103 | scDFM: Distributional Flow Matching | ICLR 2026, MMD + PAD-Transformer |
| 2602.04901 | scBIG: Module-Inductive Representations | Gene program structure for perturbation |
| 2602.10156 | STRAND: Sequence-Conditioned Transport | Genomic locus encoding for zero-shot |
| 2603.07233 | PT-RAG: Retrieval-Augmented Perturbation | Differentiable cell-type-aware retrieval |
| 2603.17380 | SCALE: Atlas-Level Endpoint Transport | BioNeMo infrastructure, 12x speedup |
| 2601.11827 | SP-FM: Shortest-Path Flow Matching | Mixture-conditioned bases for OOD |
| 2504.08328 | Conditional Monge Gap | Cross-task OT generalization |
| 2508.08312 | CFM-GP: Cross-Cell-Type Flow Matching | Single model across cell types |
| 2509.25346 | SynthPert: Synthetic Reasoning Distillation | 87% accuracy on unseen cell type |
| 2510.02903 | Cell-MNN: ODE Latent Dynamics | ICR 2026, interpretable gene interactions |
| 2410.05292 | CaLMFlow: Causal LM Flow Matching | LLM-driven flow matching |
| 2511.04790 | Uhler & Zhang: Causal Structure Survey | ICM 2026, three-problem framework |
| 2406.14302 | IEM: Identifiable Exchangeable Mechanisms | Schölkopf group, ICLR 2025 |
| 2402.09236 | Causal to Concept-Based Rep. Learning | Schölkopf group, NeurIPS 2024 |
| 2411.00614 | W1 Neural OT Solver | 25-45x speedup, ISMB/ECCB 2025 |
| 2408.10609 | PerturBench | Mode collapse finding, no dominant architecture |
| 2604.14838 | Layer-wise Foundation Model Analysis | Intermediate layers outperform final |
| 2407.10785 | SAE for Pathology Foundation Models | Monosemantic biological concepts |
| 2510.22298 | MetaCaDI: Meta-Learning Causal Discovery | Few-shot causal graph recovery |
| 2507.04704 | SPATIA: Spatial Cell Phenotypes | Multi-modal flow matching, spatial |
