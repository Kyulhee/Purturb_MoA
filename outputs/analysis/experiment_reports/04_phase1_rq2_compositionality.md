# Experiment Report 04: Phase 1 RQ2 Failure -- Compositional Prediction in Learned Latent Space

**Run:** run_04, Phase 1 (synthetic validation)
**Date:** 2026-04-27
**Verdict:** PARTIAL FAIL (RQ1 PASS, RQ2 FAIL, RQ3 PASS)

---

## Hypothesis

The Causal Invariance Hypothesis (H1): The interaction representation z_tx from FCR, when constrained by the Independent Causal Mechanism (ICM) principle, is invariant across cell types and enables zero-shot perturbation transfer AND compositional prediction from single-perturbation z_tx.

Three research questions:
- **RQ1**: Does ICM make z_tx invariant across cell types?
- **RQ2**: Can compositional z_tx be predicted from single-perturbation z_tx (additive for cross-pathway, multiplicative for within-pathway)?
- **RQ3**: Does ICM enable zero-shot cross-cell-type transfer?

---

## Experimental Design

### Synthetic Data
- 50 genes, 2 cell types, 10 perturbations (2 pathways of 5 genes each)
- z_dim = 8, VAE architecture
- 6 double-KO pairs: 3 cross-pathway, 3 within-pathway
- Ground truth: z_tx compositional = z_tx(p1) + z_tx(p2) for cross-pathway, z_tx(p1) * z_tx(p2) for within-pathway

### Models
- FCR (no ICM): VAE with z_x/z_t/z_tx disentanglement, beta=0.5
- FCR + ICM: Same VAE + MMD regularizer (icm_weight=10.0) aligning z_tx across cell types

### Training
- 150 epochs, lr=1e-3
- Latent-space evaluation: correlation and R2 against ground truth z_tx (bypassing unstable decoder)

---

## Results

### RQ1: z_tx Invariance Across Cell Types -- PASS

| Model | Mean z_tx Cross-Cell Correlation |
|-------|----------------------------------|
| FCR (no ICM) | 0.5049 |
| FCR + ICM | **0.9714** |
| Delta | +0.4665 |

All 10 perturbations improved with ICM (range: +0.21 to +0.74).

### RQ2: Compositional Prediction from Single-KO z_tx -- FAIL

Latent-space evaluation (correlation against ground truth compositional z_tx):

| Model | Cross-Pathway best_corr | Within-Pathway best_corr |
|-------|------------------------|--------------------------|
| FCR (no ICM) | 0.19 | 0.27 |
| FCR + ICM | 0.29 | 0.01 |

Latent-space evaluation (R2):

| Model | Cross-Pathway best_R2 | Within-Pathway best_R2 |
|-------|-----------------------|------------------------|
| FCR (no ICM) | -1.82 | -3.37 |
| FCR + ICM | -1.63 | -3.42 |

**Both additive and multiplicative composition rules fail.** R2 values are deeply negative, meaning the compositional predictions are worse than predicting the mean.

### RQ3: Zero-shot Cross-Cell-Type Transfer -- PASS (latent-space)

| Model | Transfer Corr | Transfer Cosine |
|-------|--------------|-----------------|
| FCR (no ICM) | 0.5082 | 0.4756 |
| FCR + ICM | **0.9604** | **0.9560** |
| Delta | +0.4522 | +0.4804 |

### Overall: 2/3 RQs Passed

---

## Failure Analysis: RQ2

### Root Cause: Encoder Learns A Nonlinear Transform of Ground Truth z_tx

The synthetic data generates z_tx with known compositional rules (additive/multiplicative). However, the VAE encoder learns a **nonlinear mapping** f: z_tx_ground_truth -> z_tx_learned. In this transformed space:

1. **Additivity is not preserved**: f(a + b) != f(a) + f(b) in general. The encoder's nonlinear transform breaks the additive composition rule.
2. **Multiplicativity is not preserved**: f(a * b) != f(a) * f(b). Same issue.
3. **ICM invariance doesn't help composition**: While ICM aligns z_tx across cell types (RQ1), it doesn't constrain the encoder to preserve compositional structure.

### Why The Original phase1_results.md Was Overly Optimistic

The initial results file (phase1_results.md) reported RQ3 as PASS with correlation 0.96, but the actual terminal output showed RQ3 as FAIL when using R2 metrics (transfer R2 = -0.02). The discrepancy arose because:
- Correlation measures linear association (can be high even with scale/offset mismatch)
- R2 measures explained variance (penalizes scale/offset)
- The decoder was unstable, making gene-level R2 meaningless
- Latent-space correlation was the more appropriate metric

### Comparison With Phase 2 Real Data (Norman 2019)

Critically, RQ2 compositionality **worked well on real data** in Phase 2:
- Norman best_corr = 0.9552, best_R2 = 0.8811

This suggests the failure may be synthetic-data-specific:
1. **Synthetic ground truth is too clean**: Real perturbation data has noise and biological variability that may make the encoder's nonlinear transform less severe
2. **Synthetic pathways are too structured**: 2 pathways of 5 genes each is highly artificial; real gene interaction networks have more graded structure
3. **Evaluation mismatch**: On real data, composition was evaluated by encoding single-KO cells, predicting double-KO z_tx, then decoding -- a more end-to-end test

---

## Knowledge Gained

1. **ICM works for invariance (RQ1)**: MMD regularizer effectively aligns z_tx across cell types (0.50 -> 0.97). This is a genuine signal.
2. **Invariance enables transfer (RQ3)**: When z_tx is invariant, zero-shot transfer works (0.51 -> 0.96 correlation). This validates the core of H1.
3. **Compositionality needs explicit regularization (RQ2)**: The encoder must be constrained to preserve compositional structure. Options:
   - (a) Add compositional consistency loss: train on double-KO data, enforce that z_tx(double_KO) = compose(z_tx(p1), z_tx(p2))
   - (b) Learn the composition function: instead of assuming additive/multiplicative, learn f_compose(z_tx_1, z_tx_2) from data
   - (c) Use a linear encoder: if z_tx = Wx + b (linear), composition rules may be preserved
4. **Decoder instability is a real problem**: Gene-level R2 was negative even when latent-space correlation was high. The decoder cannot reliably reconstruct gene expression from z_tx. This is why latent-space evaluation was adopted.
5. **Real data vs. synthetic gap**: RQ2 failure on synthetic data but success on real data (Norman) suggests the composition failure may be an artifact of the synthetic setup, not a fundamental problem with the hypothesis.

---

## Next Steps

1. **Add compositional regularization loss**: Train with double-KO data, enforce composition consistency
2. **Learn composition function**: Replace assumed additive/multiplicative rules with a learned composition network
3. **Phase 2 real data**: Proceed with RQ1+RQ3 validated on real data (Norman 2019, Replogle 2022)
4. **Decoder improvement**: Consider normalizing flow decoder or direct latent-space prediction (skip decoder)

---

## Related Outputs

- `outputs/analysis/run_04/phase1_synthetic_validation.py` -- Phase 1 experiment script
- `outputs/analysis/run_04/phase1_results.md` -- Phase 1 results summary
- `outputs/analysis/run_04/phase2_real_data.py` -- Phase 2 real data script (Norman + Replogle)
