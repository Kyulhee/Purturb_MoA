# Run 05 Results: Ablation Experiments + RQ2 Gap Analysis

**Date:** 2026-04-27 | **Status:** Completed

---

## Part A: Ablation Matrix

6 configurations tested on synthetic data (n_genes=50, 2 cell types, 10 perturbations, z_dim=8):

| # | Config | RQ1 (inv) | RQ2-cross | RQ2-within | RQ3 (transfer) |
|---|--------|-----------|-----------|------------|----------------|
| 1 | FCR baseline (no ICM) | 0.5045 | 0.1988 | 0.2457 | 0.5041 |
| 2 | FCR + ICM | 0.9548 | 0.2021 | 0.0202 | 0.9344 |
| 3 | FCR + linear z_tx head | 0.8015 | 0.4036 | 0.2481 | 0.8090 |
| 4 | FCR + ICM + linear z_tx head | 0.9660 | 0.0632 | 0.0963 | 0.9619 |
| 5 | FCR + ICM + comp consistency loss | 0.9518 | **0.7915** | 0.3050 | 0.9695 |
| 6 | FCR + ICM + linear z_tx + comp loss | **0.9918** | 0.7712 | **0.5881** | **0.9934** |

### Key findings:
- **ICM is critical for RQ1/RQ3**: baseline 0.50 -> +ICM 0.95 (invariance), 0.50 -> 0.93 (transfer)
- **ICM alone hurts RQ2**: composition drops from 0.20 to 0.02 (within-pathway)
- **Linear z_tx head alone helps moderately**: RQ2-cross 0.20 -> 0.40
- **ICM + linear z_tx paradox**: RQ1/RQ3 excellent but RQ2 worse (0.06) — ICM constrains z_tx distribution, making linear head less effective for composition
- **Compositional consistency loss is the key**: RQ2-cross jumps from 0.20 to **0.79** (config 5)
- **Full model (config 6)**: best overall — RQ1=0.99, RQ3=0.99, RQ2-cross=0.77, RQ2-within=0.59

## Part B: Encoder Nonlinearity

Linear R2 of z_tx_learned ~ z_tx_ground_truth (higher = more linear = better composition):

| Config | Linear R2 | Residual norm | Verdict |
|--------|-----------|---------------|---------|
| 1. FCR baseline | 0.6897 | 0.7210 | Highly nonlinear |
| 2. FCR + ICM | 0.8691 | 0.2007 | Moderately nonlinear |
| 3. FCR + linear z_tx head | 0.7271 | 0.6549 | Highly nonlinear |
| 4. FCR + ICM + linear z_tx | 0.8670 | 0.1968 | Moderately nonlinear |
| 5. FCR + ICM + comp loss | **0.9119** | 0.2239 | Moderately nonlinear |
| 6. FCR + ICM + linear + comp | 0.8971 | 0.2171 | Moderately nonlinear |

### Key findings:
- **ICM makes the encoder more linear**: 0.69 -> 0.87. ICM's MMD alignment constrains z_tx to be more structured
- **Comp consistency loss further increases linearity**: 0.87 -> 0.91 (best)
- **Linear z_tx head alone doesn't help**: 0.69 -> 0.73 (only marginal improvement — the upstream MLP still introduces nonlinearity)
- Even the best config (0.91) is not fully linear — some nonlinearity is intrinsic to VAE encoding

## Part C: Synthetic-Real Gap Resolution

### The core experiment

Compose single-KO z_tx (additive for cross-pathway, multiplicative for within-pathway), then:

| Evaluation space | Mean corr | Mean R2 |
|-----------------|-----------|---------|
| **Latent space** (composed z_tx vs ground truth z_tx) | 0.4757 | 0.0474 |
| **Gene space** (decoded composed z_tx vs real double-KO expression) | **0.9408** | **0.8808** |

### Explanation

**Why RQ2 fails on synthetic data (latent-space eval):**
- Ground truth z_tx is composed linearly (a+b) or multiplicatively (a*b)
- Encoder learns nonlinear mapping f(), so f(a)+f(b) != f(a+b)
- Comparing f(a)+f(b) to ground truth a+b gives poor R2

**Why RQ2 succeeds on real data (gene-space eval):**
- There is no "ground truth z_tx" — composition is evaluated end-to-end
- Encode singles -> compose in latent space -> decode -> compare to real double-KO gene expression
- The encoder-decoder pair jointly learns: decode(f(a)+f(b)) ~ decode(f(a+b))
- The decoder compensates for the encoder's nonlinearity

**Implication:** Compositionality should be evaluated in OUTPUT space (gene expression), not LATENT space. The latent space need not be compositionally structured — only the composition-after-decoding needs to match reality.

This resolves the RQ2 synthetic-real gap from run_04.

---

## Overall Conclusions

1. **ICM is essential** for invariance (RQ1) and transfer (RQ3) — confirmed across all ablation configs
2. **Compositional consistency loss is essential** for RQ2 — the only component that substantially improves composition in latent space
3. **Linear z_tx head has mixed effects** — helps alone but hurts when combined with ICM (ICM constrains distribution)
4. **The RQ2 gap is resolved**: composition should be evaluated in gene space, not latent space
   - Latent-space R2 = 0.05 vs gene-space R2 = 0.88
   - Decoder compensates for encoder nonlinearity
5. **Full model (config 6) is best overall** but config 5 (without linear head) is nearly as good for RQ2-cross and simpler

## Recommendations

1. **Evaluate composition in gene space** for all future experiments
2. **Use comp consistency loss** as standard component (configs 5 or 6)
3. **Linear z_tx head is optional** — marginal benefit, adds complexity
4. **Run configs 5-6 on Norman real data** to validate ablation on real data
