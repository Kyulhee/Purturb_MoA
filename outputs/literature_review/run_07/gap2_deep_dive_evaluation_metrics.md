# Gap 2 Deep Dive: Evaluation Metrics Don't Capture What Matters

**Date:** 2026-04-30 | **Direction:** Primary new project candidate (user-selected)

---

## 1. The Core Problem

Perturbation prediction models are optimized and compared using metrics (MSE, R², Pearson) that measure **reconstruction accuracy**, not **biological fidelity**. This creates a systematic disconnect:

- Models that score well on MSE may predict biologically wrong answers
- Models that score poorly on MSE may capture the biology correctly
- The field cannot distinguish between these cases

This is not a minor methodological inconvenience — it is a **foundational crisis** that calls into question whether the entire perturbation prediction enterprise is making progress.

---

## 2. Smoking Gun Evidence

### 2.1 Ahlmann-Eltze et al. (2025) Nature Methods — The Baseline Crisis

**Title:** "Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines"

This is the single most important paper for Gap 2. Key findings:

1. **No deep learning model outperforms an additive baseline for double perturbation prediction**
   - Tested: scGPT, scFoundation, GEARS, CPA, and others
   - Baseline: simple linear model Y = GW^T P + b
   - Result: Linear baseline wins across all benchmarks

2. **No deep learning model outperforms mean/linear baseline for unseen single perturbation prediction**
   - Even with foundation model embeddings, simple models are competitive

3. **The only advantage of DL appears with pretrained perturbation embeddings from related data**
   - Linear model + scGPT perturbation embeddings outperforms everything
   - This suggests the bottleneck is representation, not model complexity

4. **All models predominantly predict buffering interactions, rarely finding synergistic ones correctly**
   - The hardest biological predictions (synergistic GI) are exactly where models fail
   - MSE-optimal predictions converge to the mean → buffering is the "safe" prediction

5. **7+ concurrent benchmarks confirm these findings**
   - This is not a dataset artifact — it's systematic

**Critical implication:** If simple baselines outperform complex models on current metrics, there are two possible explanations:
- (a) The metrics are wrong — they reward what simple models do well (predict the mean) and penalize what complex models might do better (predict genuine biological variation)
- (b) The models are genuinely no better — deep learning adds nothing to perturbation prediction

Resolving (a) vs (b) is EXACTLY what a new evaluation framework must do.

### 2.2 SCALE / Chen et al. (2026) — The Mean-Effect Trap

**Title:** "SCALE: A foundation model for single-cell perturbation prediction via conditional transport"

Key evaluation insights:

1. **MSE incentivizes the "mean-effect trap"**
   - Models that minimize MSE tend to predict safe average profiles
   - Heterogeneous signals representing true biological responses are smoothed out
   - Quote: "evaluation protocols that overemphasize reconstruction-like accuracy while underestimating biological fidelity"

2. **Proposed alternatives: PDCorr + DE Overlap**
   - PDCorr: perturbation-directional correlation — does the predicted response vector point in the right direction?
   - DE Overlap: do the predicted differentially expressed genes overlap with observed?
   - These capture biological correctness that MSE misses

3. **Cell-Eval framework limitations**
   - Even biologically-grounded metrics (PDCorr, DE Precision, LFC Spearman, Direction Agreement) are sensitive to implementation details
   - Threshold choices (e.g., what counts as "differentially expressed") can flip conclusions

4. **Quote:** "the next bottleneck may lie less in model scaling alone and more in rebuilding the evaluation and comparison framework"

### 2.3 Shesha / Raju (2026) — Magnitude ≠ Stability

**Title:** "Geometric stability of perturbation responses reveals regulatory architecture"

Key insights:

1. **Magnitude and stability can decouple**
   - Correlation ρ=0.75-0.97 across 5 datasets — usually high
   - But discordant cases are where the biology is most interesting
   - A perturbation can produce a large but geometrically incoherent shift

2. **Pleiotropic regulators pay a "geometric tax"**
   - CEBPA, GATA1: large magnitude shifts but low directional coherence
   - These master regulators affect many pathways → response vectors scatter
   - MSE would rate these predictions poorly, but the biology (pleiotropy) is real

3. **Lineage-specific factors produce coherent shifts**
   - KLF1: large AND coherent — perturbation response is geometrically stable
   - These are the "easy" cases for any model

4. **Geometric instability independently predicts biology**
   - Associated with chaperone activation (HSPA5/BiP)
   - Instability is not noise — it reflects endoplasmic reticulum stress response
   - A metric that penalizes instability would miss this signal

5. **Metric persists in nonlinear embeddings**
   - Geometric stability measured in scGPT embedding space recapitulates gene-space findings
   - This suggests the phenomenon is embedding-invariant

---

## 3. Synthesis: Why This Gap Matters Now

### 3.1 The Convergence of Evidence

Three independent groups, using different methods, arrive at the same conclusion:

| Paper | Key Finding | Implication |
|-------|------------|-------------|
| Ahlmann-Eltze (2025) | DL models ≤ linear baselines | Current metrics can't distinguish model quality |
| SCALE (2026) | MSE creates mean-effect trap | Optimizing MSE is actively harmful |
| Shesha (2026) | Magnitude ≠ stability | Single-number summaries miss biology |

### 3.2 The Causal Chain

```
Wrong metrics → Wrong optimization target → Wrong model behavior → Wrong conclusions
     ↓                ↓                        ↓                      ↓
  MSE/R²          Mean-effect trap        Smoothing heterogeneity   "DL ≤ baseline"
                                                          ↓
                                                  Biological fidelity lost
```

The causal chain runs both ways:
- If we fix metrics → models will optimize for the right thing → DL might actually outperform baselines
- If we DON'T fix metrics → we can't tell whether DL is genuinely worse or just penalized by wrong metrics

### 3.3 Why Previous Attempts Are Insufficient

| Existing Work | What It Does | What It Misses |
|---------------|-------------|----------------|
| DE Overlap (CPA) | Binary gene set overlap | Lossy; threshold-dependent; no direction |
| PDS (ARC) | Distributional similarity | Distance metric choice dominates |
| PDCorr (SCALE) | Directional correlation | Still single-number; no gene-level resolution |
| Cell-Eval (SCALE) | Multi-metric framework | Implementation-sensitive; no theoretical grounding |
| Shesha stability | Magnitude-stability decomposition | No unified framework; single metric |
| scGraph (2026) | Graph structure preservation | Only topology; no perturbation-specific |

**No existing framework simultaneously provides:**
1. Gene-level resolution (not just perturbation-level summaries)
2. Directional coherence (not just magnitude)
3. Biological grounding (correlation with experimental outcomes)
4. Theoretical guarantees (what does the metric actually measure?)
5. Robustness to implementation choices (thresholds, embeddings, etc.)

---

## 4. Proposed Project: **BioEval — Biologically-Grounded Evaluation for Perturbation Prediction**

### 4.1 Core Research Question

**"Do current perturbation prediction metrics measure biological fidelity, and can we construct evaluation metrics that do?"**

### 4.2 Sub-Questions

**RQ1 (Diagnosis):** Which current metrics correlate with biological utility, and which are artifacts of distribution matching?
- Hypothesis: MSE/R² correlate poorly with DEG recovery, directional accuracy, and perturbation coherence
- Test: Compute full metric suite on existing models + baselines; measure metric-metric and metric-biology correlations

**RQ2 (Metric Design):** Can we construct perturbation evaluation metrics that are (a) biologically grounded, (b) gene-level resolved, (c) directionally aware, and (d) implementation-robust?
- Hypothesis: A composite metric combining directional coherence (Shesha-like), DE precision/recall, and effect-size calibration will outperform MSE/R² at predicting downstream biological utility
- Test: Propose metric; validate against expert curation + experimental outcomes

**RQ3 (Resolution of the Baseline Crisis):** Does the Ahlmann-Eltze finding (DL ≤ baseline) persist under biologically-grounded evaluation, or is it an artifact of MSE-based evaluation?
- Hypothesis: Under biologically-grounded metrics, DL models will show advantages in specific regimes (e.g., synergistic GI, cross-cell-type, high-pleiotropy perturbations) even if they underperform on MSE
- Test: Re-evaluate all models from Ahlmann-Eltze benchmark under new metrics

**RQ4 (Downstream Impact):** Do models selected by biologically-grounded metrics produce better downstream experimental outcomes (active learning, hit prioritization, cross-cell-type transfer)?
- Hypothesis: Model rankings under BioEval metrics predict downstream experimental value better than MSE/R² rankings
- Test: Compare model selection under different metrics → measure AL efficiency, precision in hit prioritization

### 4.3 Project Deliverables

1. **BioEval metric suite** — open-source Python package
   - Per-gene directional accuracy (signed error, not just magnitude)
   - DEG precision/recall at multiple thresholds (curves, not single points)
   - Geometric stability score (Shesha-like, generalized)
   - Effect-size calibration (are predicted fold-changes in the right range?)
   - Perturbation coherence (do predicted responses form biologically coherent programs?)

2. **Diagnostic re-evaluation of existing models**
   - Re-evaluate GEARS, CPA, scGPT, CPA, simple baselines under BioEval
   - Direct comparison with Ahlmann-Eltze benchmark (same models, different metrics)
   - Identify regimes where DL genuinely helps vs. where it doesn't

3. **Metric-biology correlation analysis**
   - Which metrics predict downstream utility?
   - Which metrics are redundant?
   - Which metrics reveal distinct aspects of prediction quality?

4. **Recommendations for the field**
   - Minimum reporting standard for perturbation prediction papers
   - Which metrics to report, how to compute them, what thresholds to use

### 4.4 Datasets

| Dataset | Use | Key Property |
|---------|-----|-------------|
| Replogle 2022 (K562+RPE1) | RQ1-3: Primary evaluation | 848 shared perturbations, cross-cell-type |
| Norman 2019 | RQ1-3: Combinatorial evaluation | 128 double-KOs with ground-truth GI |
| Ahlmann-Eltze benchmark suite | RQ3: Direct replication | 7+ benchmarks with existing model outputs |
| PORTAL 2026 | RQ3-4: Large-scale validation | 665K pairwise perturbations |
| Tahoe-100M | RQ4: Benchmark platform | 100M cells, diverse perturbations |

### 4.5 Baselines and Comparisons

| Method | Why Include |
|--------|------------|
| Additive baseline (linear) | Ahlmann-Eltze winner — does it still win under BioEval? |
| GEARS | GNN-based, most cited combinatorial predictor |
| CPA | Compositional autoencoder, standard comparison |
| scGPT + linear | Best performer in Ahlmann-Eltze with pretrained embeddings |
| SCALE | Foundation model, proposed PDCorr |
| Mean predictor | Trivial baseline — should score poorly on any reasonable metric |

### 4.6 Key Metrics to Compare

| Metric Category | Specific Metrics | What It Captures |
|-----------------|-----------------|------------------|
| Distribution matching | MSE, MAE, R², Pearson | Current standard |
| DEG recovery | DE overlap, DE precision, DE recall, F1 | Gene-level correctness |
| Directional | PDCorr, direction agreement, Shesha stability | Coherence of response |
| Effect-size | LFC Spearman, fold-change correlation | Magnitude accuracy |
| Composite | BioEval (proposed) | Unified biological fidelity |

---

## 5. Novelty Assessment

### 5.1 What Exists

- **Cell-Eval (SCALE, 2026)**: Multi-metric framework, but implementation-sensitive and no theoretical analysis of metric properties
- **Shesha (2026)**: Single new metric (geometric stability), no unified framework
- **Ahlmann-Eltze (2025)**: Diagnostic (metrics are broken), no proposed fix
- **ARC Virtual Cell Challenge**: PDS as primary metric, limited scope

### 5.2 What's New in This Project

| Contribution | Direct Competitors | Novelty Level |
|-------------|-------------------|---------------|
| Unified biologically-grounded evaluation suite | Cell-Eval (partial) | High — no existing suite has all 5 properties |
| Diagnostic: does metric choice flip model rankings? | Ahlmann-Eltze (diagnosis only) | High — no one has shown this quantitatively |
| Resolution of baseline crisis under new metrics | No one | Very High — the field's most pressing open question |
| Metric-biology correlation (which metrics predict utility?) | No one | Very High — no prior work |
| Minimum reporting standard for the field | No one | Moderate — community resource |

### 5.3 Risk: This Could Be a Negative Result

If DL models STILL underperform baselines under biologically-grounded metrics, that's actually an important finding — it means the baseline crisis is real, not a metric artifact. This is still publishable (Nature Methods level, following Ahlmann-Eltze).

If DL models outperform baselines under new metrics, it resolves the crisis and provides the field with better tools. Also high-impact.

**Either outcome is publishable.** This is a rare low-risk/high-reward project.

---

## 6. Comparison with Gap 1 (Conformal Prediction)

| Dimension | Gap 1 (Conformal) | Gap 2 (Evaluation) |
|-----------|-------------------|---------------------|
| Novelty | High (no prior work) | High (no unified framework) |
| Risk | Medium (exchangeability may not hold) | Low (any outcome is publishable) |
| Impact breadth | UQ users only | Entire perturbation prediction field |
| Feasibility | High | High |
| User interest | — | **User selected this direction** |
| Dependency | Requires base model to calibrate | Model-agnostic |
| Timeline | Medium | Short-Medium |
| Synergy with prior work | Moderate (we found coverage problem) | High (our own CPA>FCR finding is evidence for this gap) |

---

## 7. Concrete First Steps

1. **Week 1:** Re-implement Ahlmann-Eltze benchmark suite. Compute all models' predictions. Compute MSE/R² baseline rankings.

2. **Week 2:** Implement BioEval metric suite (v0.1). Compute on all model predictions. Compare rankings.

3. **Week 3:** Metric-biology correlation analysis. Which metrics predict DEG recovery, directional accuracy, and perturbation coherence?

4. **Week 4:** Write up results. If rankings flip under BioEval → high-impact paper. If rankings don't flip → important negative result paper.

---

## 8. Key Papers

| Paper | arXiv/DOI | Role |
|-------|-----------|------|
| Ahlmann-Eltze et al. (2025) | 10.1038/s41592-025-02772-6 | Smoking gun: DL ≤ baseline |
| Chen et al. (2026) SCALE | arXiv:2603.17380 | Mean-effect trap, PDCorr, Cell-Eval |
| Raju (2026) Shesha | arXiv:2604.16642 | Geometric stability, magnitude≠direction |
| Roohani et al. (2023) GEARS | Nature Biotech | Standard model to re-evaluate |
| Lotfollahi et al. (2023) CPA | Mol Syst Biol | Standard model to re-evaluate |
| Tang & Norman (2026) PORTAL | bioRxiv | Large-scale validation resource |
