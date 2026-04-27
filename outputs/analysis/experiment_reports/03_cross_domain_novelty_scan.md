# Experiment Report 03: Cross-Domain Novelty Scan -- Hypothesis Selection Process

**Run:** run_04 (early phase)
**Date:** 2026-04-27
**Verdict:** TRANSITION (not a failed experiment, but a filtered selection process)

---

## Context

After run_01-02 showed that GNN+XGBoost and AL for metabolic modeling lacked novelty, the research direction pivoted to Perturb-seq / single-cell perturbation prediction. A systematic cross-domain novelty scan was conducted to identify high-novelty hypothesis candidates.

---

## Scan Process

### Method
- Systematic arXiv search (2024-2026) for 10 cross-domain intersections with perturbation prediction
- Key paper deep-reads for each domain
- Gap analysis: what exists vs. what is novel/untested
- Assessment on novelty, feasibility, and impact

### 10 Domains Scanned

| # | Domain | Novelty | Feasibility | Impact | Verdict |
|---|--------|---------|-------------|--------|---------|
| 1 | Causal Representation Learning | HIGH | Medium | Very High | SELECTED (became H1) |
| 2 | In-Context Learning / Foundation Models | MEDIUM | High | Very High | Rejected -- MapPFN already implements core idea |
| 3 | Meta-Learning / Few-Shot | HIGH | Medium | High | Considered -- folded into H2 |
| 4 | Equivariant Neural Networks | HIGH | Low-Medium | Medium | Rejected -- feasibility too low, correct symmetry group unclear |
| 5 | Neural ODE / CNF | MEDIUM | Medium | High | Rejected -- flow matching already well-explored |
| 6 | Mechanistic Interpretability | HIGH | Medium | High | Selected as H3 (audit hypothesis) |
| 7 | Optimal Transport | LOW | High | Medium | Rejected -- saturated domain (CellOT, W1 solver, SCALE) |
| 8 | Compositional Generalization | HIGH | Medium | Very High | SELECTED (became core of H1) |
| 9 | Physics-Informed NNs | HIGH | Medium | High | Deferred -- requires high-quality GRN databases |
| 10 | IRM / OOD Generalization | MEDIUM | Medium | High | SELECTED (merged into H1 as ICM invariance) |

---

## Top 3 Hypotheses Selected

### H1: Causal Invariance Hypothesis (Domains 1 + 8 + 10)
"The interaction representation z_tx from FCR, when constrained by the Independent Causal Mechanism principle, is invariant across cell types and enables zero-shot perturbation transfer."

- Combines: FCR disentanglement + ICM invariance + compositional prediction
- Why selected: Highest novelty/feasibility/impact combination; FCR provides the representation, ICM provides the invariance objective, composition provides the application

### H2: Compositional Module Hypothesis (Domains 3 + 8)
"Gene perturbation effects decompose into pathway-level modules that compose via a learnable interaction function."

- Merged into H1's composition component
- Not pursued separately

### H3: Mechanistic Interpretability Audit (Domain 6)
"Current perturbation prediction models learn statistical shortcuts, not biological mechanisms."

- Deferred -- requires training multiple baseline models first
- Could be a future validation step for H1

---

## Why 7 Domains Were Rejected

| Domain | Rejection Reason |
|--------|-----------------|
| In-Context Learning | MapPFN (Jan 2026) already implements core idea; novelty = MEDIUM |
| Equivariant NNs | Feasibility LOW -- correct biological symmetry group is unknown; partial-order equivariance is theoretically appealing but hard to implement |
| Neural ODE / CNF | Flow matching is already well-explored for perturbation (scDFM, CFM-GP, SP-FM); the specific ODE formulation would need time-series data which Norman 2019 lacks |
| Optimal Transport | Saturated -- CellOT, W1 solver, Conditional Monge Gap, SCALE all exist; incremental improvements only |
| Physics-Informed NNs | Requires high-quality pathway/GRN databases; KEGG/Reactome are incomplete; soft constraints may be more practical but less principled |

---

## Knowledge Gained

1. **Novelty requires timing**: MapPFN (Jan 2026) pre-empted the in-context learning hypothesis. Rapid literature review is essential before committing to a direction.
2. **Feasibility filters matter**: Equivariant NNs scored HIGH novelty but LOW feasibility. Without a clear symmetry group, implementation would be speculative.
3. **Domain intersection is powerful**: H1 emerged from combining 3 domains (causal representations + compositionality + IRM). None alone was as compelling as the combination.
4. **Perturb-seq is the right domain**: 20+ methods published in 2025-2026, yet fundamental problems remain (compositional prediction, cross-cell-type transfer). This is an active frontier where high-impact work is possible.

---

## Related Outputs

- `outputs/analysis/run_04/cross_domain_novelty_scan.md` -- Full 10-domain scan with detailed analysis
- `outputs/analysis/run_04/perturb_seq_unsolved_problems_review.md` -- 5 core unsolved problems synthesis
