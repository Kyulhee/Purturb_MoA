# Phase 1: Synthetic FCR+ICM Validation Results

## Configuration
- 50 genes, 2 cell types, 10 perturbations (2 pathways), z_dim=8
- 150 epochs, lr=1e-3, beta=0.5, icm_weight=10.0
- Double-KO pairs: 3 cross-pathway, 3 within-pathway

## Results

### RQ1: z_tx Invariance Across Cell Types — PASS
- FCR (no ICM): mean z_tx cross-cell correlation = 0.5049
- FCR + ICM: mean z_tx cross-cell correlation = 0.9714
- Delta: +0.4665
- All 10 perturbations improved with ICM (range: +0.21 to +0.74)

### RQ2: Compositional Prediction from Single-KO z_tx — FAIL
Latent-space evaluation (correlation against ground truth compositional z_tx):
- FCR (no ICM): cross-pathway best_corr=0.19, within-pathway best_corr=0.27
- FCR + ICM: cross-pathway best_corr=0.29, within-pathway best_corr=0.01
- Composition rules (additive for cross-path, multiplicative for within) don't hold in learned z_tx space
- Possible cause: encoder learns nonlinear transform of ground truth z_tx

### RQ3: Zero-shot Cross-Cell-Type Transfer — PASS
Latent-space evaluation (correlation between source and target z_tx):
- FCR (no ICM): transfer corr=0.5082, cos=0.4756
- FCR + ICM: transfer corr=0.9604, cos=0.9560
- Transfer improvement: +0.4522
- ICM makes z_tx nearly invariant, enabling near-perfect cross-cell-type transfer

## Overall: 2/3 RQs passed — SIGNAL DETECTED

## Key Insights
1. **ICM works for invariance**: MMD regularizer effectively aligns z_tx across cell types (0.50 -> 0.97)
2. **Invariance enables transfer**: When z_tx is invariant, zero-shot transfer works (0.51 -> 0.96)
3. **Compositionality needs work**: Assumed additive/multiplicative rules don't apply in learned latent space
   - The encoder may learn a rotated/transformed representation where composition isn't simple
   - Need to either: (a) learn the composition function, or (b) add compositional regularization during training

## Next Steps
- RQ2 fix: Add compositional consistency loss (train on double-KO data, enforce composition)
- Phase 2: Test on real data (Norman 2019, Replogle 2022) with RQ1+RQ3 validated
