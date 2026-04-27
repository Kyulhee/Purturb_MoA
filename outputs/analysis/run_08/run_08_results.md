# Run 08 Results: Baseline Comparison for Cross-Cell-Type Transfer

**Date:** 2026-04-27 | **Status:** Completed (GEARS failed)

---

## Summary

Compared FCR-ICM against baselines on RQ3 (zero-shot K562→RPE1 transfer, 843 shared perturbations).

| Method | R2 | Corr | N |
|--------|-----|------|---|
| Mean Shift | 0.8189 | 0.9125 | 843 |
| FCR (no ICM) | -0.1690 | 0.4262 | 843 |
| **FCR + ICM (ours)** | **0.9210** | **0.9658** | 843 |
| GEARS | Failed (API error) | — | — |

## Key Findings

### 1. Mean Shift is a surprisingly strong baseline (R2=0.82)

The simplest possible transfer method — adding K562's perturbation delta to RPE1's control mean — achieves R2=0.82 and corr=0.91. This means perturbation effects are largely additive and cell-type-invariant at the mean level.

However, Mean Shift has critical limitations:
- Only predicts population-level means (no single-cell resolution)
- Cannot capture nonlinear/interaction effects
- Assumes perfect additivity of perturbation effects
- No uncertainty quantification

### 2. FCR+ICM significantly outperforms Mean Shift (+0.10 R2)

FCR+ICM R2=0.921 vs Mean Shift R2=0.819 (+0.102). The gap is meaningful:
- FCR+ICM captures single-cell heterogeneity (generative model)
- ICM explicitly enforces z_tx invariance rather than assuming it
- FCR+ICM can generalize to unseen perturbation combinations (RQ2)

### 3. FCR without ICM fails dramatically (R2=-0.17)

Without ICM, the model learns cell-type-specific z_tx that transfers poorly (corr=0.43). This confirms ICM is essential, not just helpful.

### 4. GEARS evaluation failed

GEARS API requires `pd_k562.prepare_split()` + `pd_k562.get_dataloader()` before `GEARS(pd_k562)` initialization. Fix:
```python
pd_k562.prepare_split(split='simulation', seed=1)
pd_k562.get_dataloader(batch_size=32)
model = GEARS(pd_k562, device='cpu')
```
GEARS does NOT natively support cross-cell-type transfer. Workaround: swap `model.ctrl_adata` with RPE1 ctrl cells.

---

## Literature Review: Cross-Cell-Type Perturbation Transfer

### Directly Competing Methods

1. **C3TL (Scholkemper & Mukherjee, 2026, arXiv:2603.13051)**
   - "Causal Cellular Context Transfer Learning"
   - **Directly competing approach**: uses causal invariance for cross-cell-type perturbation transfer
   - Lightweight framework, requires only bulk molecular data
   - Claims competitive with SOTA foundation models with much smaller model sizes
   - **Key overlap with FCR-ICM**: causal invariance principle, context transfer, lightweight architecture
   - **Key difference**: C3TL uses bulk data; FCR-ICM uses single-cell data + VAE decomposition

2. **STRAND (Fu et al. 2026, arXiv:2602.10156)**
   - "Sequence-Conditioned Transport for Single-Cell Perturbations"
   - Conditions on regulatory DNA sequence at perturbation locus
   - Tested on K562/Jurkat/RPE1 — same cell lines as our Replogle dataset
   - Improves transfer to novel cell lines by up to 0.14 in Pearson correlation
   - **Key difference**: sequence-conditioned; FCR-ICM is perturbation-embedding-conditioned

3. **Wang et al. 2026 (Nature Computational Science, DOI:10.1038/s43588-025-00887-6)**
   - "Predicting drug responses of unseen cell types through transfer learning with foundation models"
   - Uses scGPT/foundation models for cross-cell-type drug response prediction
   - **Key difference**: drug perturbation (not genetic); foundation model approach

4. **Wei et al. 2026 (Nature Methods, DOI:10.1038/s41592-025-02980-0)**
   - "Benchmarking algorithms for generalizable single-cell perturbation response prediction"
   - Comprehensive benchmark of perturbation prediction algorithms
   - Tests generalization to unseen perturbations, unseen cell types, unseen combinations
   - 13 citations — important reference for our paper's related work section

5. **dbDiffusion (Shang, Wei & Roeder, 2025, PNAS, DOI:10.1073/pnas.2525268122)**
   - Diffusion-based debiasing framework for perturbation prediction
   - Does NOT rely on foundation models
   - Transfers information across related experimental conditions
   - **Key difference**: diffusion model; no explicit causal invariance

6. **scPerb (Tang et al. 2025, J Adv Res, DOI:10.1016/j.jare.2024.10.035)**
   - Style transfer-based VAE for perturbation prediction
   - **Key overlap**: VAE-based, style transfer concept similar to ICM alignment

### Key Insight: C3TL is the Most Direct Competitor

C3TL (arXiv:2603.13051, March 2026) independently arrived at the same core idea as FCR-ICM:
- Causal invariance for cross-cell-type perturbation transfer
- Lightweight model (vs foundation models)
- Claims SOTA-competitive performance

**Differences to emphasize in our paper:**
1. FCR-ICM uses single-cell VAE with z_x/z_t/z_tx decomposition (finer-grained than C3TL)
2. FCR-ICM validates on genetic perturbations (CRISPR); C3TL may differ in perturbation type
3. FCR-ICM explicitly demonstrates the decoder compensation effect (latent vs gene space evaluation)
4. FCR-ICM shows compositional prediction (RQ2) in addition to transfer (RQ3)
5. Our Mean Shift baseline (R2=0.82) provides important context for interpreting transfer results

### Critical Implication: Mean Shift R2=0.82 Weakens the "Novelty" Argument

The strong Mean Shift performance suggests that **simple additive transfer of perturbation effects already works well at the population level**. To establish FCR+ICM's unique contribution, we need to demonstrate:

1. **Single-cell resolution**: Mean Shift only predicts means; FCR+ICM generates full distributions
2. **Compositional prediction**: Mean Shift cannot predict unseen combinations; FCR+ICM can (RQ2)
3. **Non-additive effects**: For perturbations with cell-type-specific effects, Mean Shift would fail
4. **Per-cell prediction**: Individual cell-level predictions, not just population means

---

## Additional Baseline Methods (from literature research)

### CPA (Compositional Perturbation Autoencoder)
- **Cross-cell-type: Native support** via adversarial disentanglement
- Installation: `pip install cpa-tools` (v0.8.1)
- Uses adversarial training to disentangle perturbation effect from cell-type latent
- Train on K562 (all) + RPE1 (ctrl only), mark RPE1 perturbed as OOD
- **Priority: HIGH** — most principled baseline for cross-cell-type comparison

### scGen
- **Cross-cell-type: Native support** via delta vector
- Installation: `pip install scgen` (v2.1.0, depends on scvi-tools)
- `model.predict(celltype_to_predict='RPE1')` — directly supports target cell type
- Learns a single linear delta vector per perturbation (simplest baseline)
- **Priority: MEDIUM** — simplest baseline, but linear only

### GEARS (fix needed)
- **Cross-cell-type: Workaround only** — swap ctrl_adata
- Fix: call `pd_k562.prepare_split()` + `pd_k562.get_dataloader()` before `GEARS(pd_k562)`
- GEARS has NO native cross-cell-type support (README explicitly states this)
- The workaround is equivalent to Mean Shift in principle (perturbation invariance assumption)
- **Priority: HIGH** — main SOTA comparison target

### CellOT
- Not pip-installable, script-based, no CRISPR support
- **Priority: LOW**

---

## Next Steps

1. **Fix GEARS baseline**: Call `prepare_split()` + `get_dataloader()` before `GEARS()` initialization
2. **Add CPA baseline**: `pip install cpa-tools`, train with cell_type as covariate, evaluate OOD prediction
3. **Add scGen baseline**: `pip install scgen`, simplest linear delta baseline
4. **Analyze Mean Shift failure cases**: Which perturbations does Mean Shift fail on? Are they non-additive?
5. **Single-cell evaluation**: Compare Mean Shift vs FCR+ICM at single-cell level (not just means)
6. **Read C3TL in detail**: Understand exact methodology overlap and differences
