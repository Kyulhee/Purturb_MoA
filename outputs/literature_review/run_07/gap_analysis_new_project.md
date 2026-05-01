# Perturb-seq Gap Analysis for New Project

**Date:** 2026-04-30 | **Purpose:** Identify unmet needs for a new perturb-seq research project

---

## 1. Field Landscape (2024-2026)

### 1.1 Perturbation Prediction Methods — Rapidly Crowding

| Method | Year | Key Idea | Limitation |
|--------|------|----------|------------|
| GEARS | 2023 | GNN + GRN for GI classification | Requires combinatorial training data |
| CPA | 2023 | Compositional autoencoder | OOD failure (DEG 0.85→0.38) |
| CellOT | 2023 | Neural optimal transport | Cross-patient only, no gene perturbation |
| scDFM | 2026 | Distributional flow matching | No combinatorial structure |
| STRAND | 2026 | Sequence-conditioned transport | Locus-aware but no GI detection |
| SCALE | 2026 | LLaMA-based + conditional transport | Foundation model, evaluation gap noted |
| PerturbDiff | 2026 | Functional diffusion | Unpaired data challenge |
| Departures | 2025 | Schrödinger bridge | Distribution transport, no UQ |
| C3TL | 2026 | Causal context transfer | Competes with FMs on efficiency |
| PRESCRIBE | 2025 | Bayesian evidential regression | **Only UQ method for perturbation** |
| scPPDM | 2025 | Diffusion for drug response | Drug-specific, not gene perturbation |

### 1.2 Key Observations
- **20+ methods in 2 years** — the field is crowded for point prediction
- **Foundation models arriving**: SCALE (LLaMA-based), scGPT embeddings used in Shesha
- **ARC Virtual Cell Challenge** driving benchmark standardization
- **PORTAL (2026)**: 665,856 pairwise perturbations, 46M clones — massive validation resource now available
- **Tahoe-100M**: New 100M-cell benchmark dataset

---

## 2. Identified Gaps (Ranked by Need × Feasibility)

### Gap 1: Uncertainty Calibration is Broken and Unsolved ⭐⭐⭐⭐⭐

**The Problem:**
- Every perturbation prediction model outputs point estimates
- PRESCRIBE (2025) is the **only** method addressing UQ — uses evidential deep learning
- But PRESCRIBE provides only a "confidence score" — NOT calibrated prediction intervals
- Our own run_12 showed: MC Dropout 90% CI covers only **29.3%** of genes (target: 85-95%)
- Norman dataset: coverage = 23.5%
- This is **severe under-coverage** — model uncertainty is miscalibrated

**Why It Matters:**
- Biologists using these models need to know WHICH predictions to trust
- A perturbation prediction with no error bar is scientifically incomplete
- Without calibration, you cannot do risk-aware experiment design (active learning is on shaky ground)
- Conformal prediction exists as a general framework but has **never been applied to single-cell perturbation prediction**

**Evidence of Need:**
- PRESCRIBE paper itself shows: filtering untrustworthy predictions gives 3%+ improvement — UQ has practical value
- Our run_12: rho_mc=0.786 (MC Dropout correlates well with error) but coverage=0.293 (intervals are wrong size)
- The correlation is there; the calibration is what's missing

**Direct Competitors:** 0 (no one has done conformal prediction for perturbation response)

**Feasibility:** High
- Conformal prediction is distribution-free — works with any base model
- Holdout calibration sets are natural (held-out perturbations)
- Exchangeability assumption: perturbations are approximately exchangeable
- Can build on top of existing models (GEARS, CPA, etc.)

---

### Gap 2: Evaluation Metrics Don't Capture What Matters ⭐⭐⭐⭐

**The Problem:**
- SCALE (2026) explicitly states: "evaluation protocols that overemphasize reconstruction-like accuracy while underestimating biological fidelity"
- Standard metrics: MSE, R², Pearson — measure distribution matching, not biological correctness
- DEG overlap is binary and lossy
- Shesha (2026, Raju): introduces geometric stability (directional coherence) as complementary to effect magnitude — magnitude and stability can decouple
- scGraph (2026): biological structure preservation

**Why It Matters:**
- Models optimize for MSE but biologists care about: "Did you get the right genes changing? In the right direction? Is the perturbation coherent?"
- A model with R²=0.95 but wrong DEG ranking is useless for hit prioritization
- The "simple baseline" problem (CPA > FCR in our run_12) may be an artifact of wrong metrics

**Evidence of Need:**
- Shesha: stability correlates ρ=0.75-0.97 with magnitude across 5 datasets — but discordant cases reveal real biology (pleiotropic regulators)
- SCALE: PDCorr + DE overlap as biologically meaningful alternatives
- Our own finding: prod formula rho=0.437 (PASS) vs A7 rho=0.326 (PARTIAL) — metric choice changes conclusion

**Feasibility:** Moderate
- Requires defining what "biologically correct" means — domain-specific
- PORTAL provides ground truth for validation at unprecedented scale
- Tahoe-100M as benchmark platform

---

### Gap 3: Out-of-Distribution Generalization is the Real Bottleneck ⭐⭐⭐⭐

**The Problem:**
- All models trained on seen perturbations/cell types
- Cross-cell-type transfer: our project showed it works partially (rho=0.33-0.44) but ICM doesn't help
- C3TL (2026): claims competitive with foundation models using causal architecture — but still limited
- CPA OOD failure: DEG 0.85→0.38
- No method provides reliable predictions for truly novel perturbation combinations

**Why It Matters:**
- The whole point of perturbation prediction is to predict what you HAVEN'T measured
- If you need training data for every combination, you might as well do the experiment
- Cross-cell-type is critical for translational applications (cell line → primary cell)

**Evidence of Need:**
- Hetzel et al. (2022): explicit OOD evaluation shows dramatic degradation
- Our run_12: cross-CT transfer works but with moderate correlation (0.33-0.44)
- C3TL: "such [foundation model] approaches are computationally expensive and may not always generalize well"

**Feasibility:** Moderate
- Domain adaptation/transfer learning is well-studied in ML
- Perturb-seq has natural OOD structure (unseen perturbations, unseen cell types)
- Conformal prediction under covariate shift exists (Guan 2021, Jonkers 2024)

---

### Gap 4: No Principled Active Learning for Experiment Design ⭐⭐⭐

**The Problem:**
- NAIAD (2024): adaptive embedding + AL for combinatorial screening, up to 40% improvement
- Our run_12: AL 2.67× improvement for epistasis detection
- But: AL requires reliable uncertainty estimates, which are miscalibrated (Gap 1)
- Current AL is heuristic — no formal optimality guarantees

**Why It Matters:**
- Perturb-seq experiments are expensive (~$5K-50K per run)
- Intelligent experiment design can dramatically reduce cost
- AL + calibrated UQ = principled "where to measure next"

**Evidence of Need:**
- NAIAD: 40% improvement in combinatorial screening efficiency
- Our project: AL works (2-5× improvement) but on uncalibrated uncertainty
- No benchmark or formal framework for perturb-seq AL exists

**Feasibility:** Moderate-High
- AL theory is mature (Bayesian optimization, bandits)
- Conformal prediction can provide calibrated acquisition functions
- Natural setup: sequential batch design

---

### Gap 5: Combinatorial Perturbation Prediction Remains Unsolved ⭐⭐⭐

**The Problem:**
- GEARS requires combinatorial training data
- CPA fails on OOD combinations
- scBIG: 6.7% improvement with module structure — incremental
- No method can reliably predict double-KO effects from single-KO data alone
- The "additivity assumption" holds for 86% of Norman double-KOs but fails for the 14% that matter most (synthetic lethals, buffering)

**Why It Matters:**
- Drug combination therapy requires combinatorial prediction
- Synthetic lethal screening (cancer) is exactly the regime where additivity fails
- PORTAL's 665K pairwise data provides ground truth at scale

**Evidence of Need:**
- Norman et al.: CBL/CNN1 synergy cannot be predicted from singles
- Our project: epistasis is partially conserved across cell types — but signal is moderate
- Valenzuela (2025): formula choice changes GI classification — the definition itself is unstable

**Feasibility:** Low-Moderate
- Deep problem — biology is inherently non-additive
- PORTAL enables validation but not necessarily solution
- May need fundamentally different modeling approach (mechanistic + data-driven)

---

## 3. Recommended New Project Direction

### Primary Recommendation: **BioEval — Biologically-Grounded Evaluation for Perturbation Prediction** ⭐ USER SELECTED

**Core idea:** Construct evaluation metrics for perturbation prediction that measure biological fidelity rather than reconstruction accuracy, and diagnose whether the "DL ≤ baseline" crisis is a metric artifact or a genuine limitation.

**Why this gap (user-selected):**
1. **Most pressing open question in the field**: Ahlmann-Eltze (2025) Nature Methods showed DL ≤ linear baseline — but is this because metrics are wrong, or models are genuinely no better?
2. **Converging evidence from 3 independent groups**: Ahlmann-Eltze (baseline crisis), SCALE (mean-effect trap), Shesha (magnitude≠stability)
3. **Low risk / high reward**: Either outcome is publishable. If rankings flip under BioEval → resolves crisis. If rankings don't flip → important negative result (Nature Methods level).
4. **Our unique experience**: Our own CPA>FCR finding (run_12: 0.430 vs 0.367) is micro-evidence for this gap. We understand the metric sensitivity problem from 3-formula sensitivity analysis.
5. **Impact breadth**: Entire perturbation prediction field needs this. Every new method paper will use the evaluation framework.
6. **Feasibility is very high**: No new model training needed. Re-evaluate existing model outputs under new metrics.

**Potential contributions:**
1. BioEval metric suite — gene-level, directional, biologically grounded, implementation-robust
2. Diagnostic: does metric choice flip model rankings? (No prior quantitative analysis)
3. Resolution of baseline crisis under biologically-grounded metrics (field's most pressing question)
4. Metric-biology correlation analysis — which metrics predict downstream utility?
5. Minimum reporting standard for perturbation prediction papers

**Detailed proposal:** See `outputs/literature_review/run_07/gap2_deep_dive_evaluation_metrics.md`

**Datasets:**
- Replogle 2022 (K562+RPE1, 848 shared perturbations) — cross-cell-type
- Norman 2019 (128 double-KOs) — combinatorial
- Ahlmann-Eltze benchmark suite (7+ benchmarks) — direct replication
- PORTAL 2026 (665K pairwise) — large-scale validation
- Tahoe-100M — benchmark platform

### Alternative Direction: **Calibrated Uncertainty for Perturbation Prediction (Gap 1)**

Apply conformal prediction to single-cell perturbation models. No direct competitors. Higher novelty but narrower impact (UQ users only). Our coverage failure (0.293) provides unique entry point.

### Alternative Direction: **Cross-Cell-Type Transfer with Evaluation Guarantees (Gap 2 + Gap 3)**

Combine evaluation reform with cross-cell-type transfer: evaluate whether models that perform well under BioEval also transfer better across cell types. More specific, higher integration potential.

---

## 4. Risk Assessment

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Conformal prediction exchangeability violation | Medium | High | Use perturbation-conditional conformal; weighted conformal under covariate shift |
| Coverage intervals too wide to be useful | Medium | Medium | Adaptive coverage levels; per-perturbation calibration |
| Base model choice affects results | Low | Medium | Test on multiple base models (CPA, GEARS, simple baseline) |
| Reviewer unfamiliarity with conformal prediction | Low | Medium | Clear tutorial-style presentation; biological motivation first |
| PORTAL/Tahoe data access issues | Low | High | Start with Replogle+Norman; PORTAL as validation |

---

## 5. Key Papers for New Project

| Paper | Relevance | Key Takeaway |
|-------|-----------|--------------|
| PRESCRIBE (Cheng et al., 2025, arXiv:2510.07964) | Direct competitor (UQ) | Evidential regression for perturbation UQ; confidence score only, not calibrated |
| SCALE (Chen et al., 2026, arXiv:2603.17380) | Evaluation gap | "Evaluation overemphasizes reconstruction"; PDCorr + DE overlap as alternatives |
| Shesha (Raju, 2026, arXiv:2604.16642) | Evaluation metric | Geometric stability as complementary to magnitude; pleiotropic regulators |
| C3TL (Scholkemper & Mukherjee, 2026, arXiv:2603.13051) | Cross-cell-type | Causal architecture for context transfer; competitive with FMs |
| Angelopoulos & Bates (2021, arXiv:2107.07511) | Conformal prediction tutorial | Distribution-free UQ framework |
| Guan (2021, arXiv:2106.08460) | Localized conformal | Per-sample adaptive calibration |
| Jonkers et al. (2024, arXiv:2404.15018) | Conformal under covariate shift | Extends CP to non-IID settings |
| Valenzuela et al. (2025) | GI formula sensitivity | Product neutrality in biological systems |
| NAIAD (Qin et al., 2024, arXiv:2411.12010) | Active learning for perturbation | AL for combinatorial screening; 40% improvement |
