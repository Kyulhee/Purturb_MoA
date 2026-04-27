"""
Phase 1: Synthetic FCR+ICM Validation
=======================================
Fast proof-of-concept: generate synthetic perturbation data with KNOWN invariant z_tx,
train FCR ± ICM, test if ICM improves cross-cell-type invariance and compositionality.

Key questions:
  RQ1: Does ICM regularization make z_tx more invariant across cell types?
  RQ2: Can we compose z_tx from single perturbations to predict double-KO effects?
  RQ3: Does invariant z_tx enable zero-shot cross-cell-type transfer?
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal
from sklearn.metrics import r2_score
import warnings
warnings.filterwarnings("ignore")

torch.manual_seed(42)
np.random.seed(42)

# ============================================================
# 1. Synthetic Data Generator (ground truth known)
# ============================================================

class SyntheticPerturbationData:
    """
    Generates single-cell perturbation data with known latent structure.

    Ground truth:
      z_x ~ N(mu_c, Sigma_c)  -- cell-type specific covariate
      z_t = embedding(perturbation_id)  -- treatment embedding
      z_tx = f_interaction(z_t, pathway_id)  -- INVARIANT across cell types

    Observation:
      x = decoder(z_x, z_t, z_tx) + noise

    Compositionality:
      - Within-pathway double KO: multiplicative (z_tx_{A+B} = z_tx_A * z_tx_B)
      - Cross-pathway double KO: additive (z_tx_{A+B} = z_tx_A + z_tx_B)
    """

    def __init__(self, n_genes=100, n_cell_types=2, n_perturbations=10,
                 n_cells_per_condition=200, z_dim=10):
        self.n_genes = n_genes
        self.n_cell_types = n_cell_types
        self.n_perturbations = n_perturbations
        self.n_cells = n_cells_per_condition
        self.z_dim = z_dim

        # Two pathways of 5 perturbations each
        self.n_pathways = 2
        self.perturbations_per_pathway = n_perturbations // self.n_pathways
        self.pathway_assignment = {}
        for i in range(n_perturbations):
            self.pathway_assignment[i] = i // self.perturbations_per_pathway

        # Cell type-specific covariate distributions
        self.cell_type_means = np.random.randn(n_cell_types, z_dim) * 2.0
        self.cell_type_stds = np.abs(np.random.randn(n_cell_types, z_dim)) * 0.5 + 0.5

        # Perturbation embeddings (z_t) -- shared across cell types
        self.z_t_ground = np.random.randn(n_perturbations, z_dim) * 0.5

        # Interaction effects (z_tx) -- INVARIANT by design
        # Each perturbation has a pathway-level module effect
        self.z_tx_ground = np.random.randn(n_perturbations, z_dim) * 1.0

        # Decoder weights (shared across cell types)
        self.W_dec = np.random.randn(3 * z_dim, n_genes) * 0.3
        self.b_dec = np.random.randn(n_genes) * 0.1

        # Combinatorial ground truth
        self.double_ko_pairs = []
        # 3 cross-pathway pairs (additive)
        for i in range(3):
            p1 = i  # pathway 0
            p2 = self.perturbations_per_pathway + i  # pathway 1
            self.double_ko_pairs.append((p1, p2))
        # 3 within-pathway pairs (multiplicative)
        for i in range(3):
            p1 = i  # pathway 0
            p2 = i + 3  # also pathway 0
            if p2 < self.perturbations_per_pathway:
                self.double_ko_pairs.append((p1, p2))
        # If we didn't get enough within-pathway, add more
        while len([p for p in self.double_ko_pairs
                   if self.pathway_assignment[p[0]] == self.pathway_assignment[p[1]]]) < 3:
            i = np.random.randint(0, self.perturbations_per_pathway)
            j = np.random.randint(0, self.perturbations_per_pathway)
            if i != j and (i, j) not in self.double_ko_pairs:
                self.double_ko_pairs.append((i, j))

    def _decoder(self, z_x, z_t, z_tx):
        """Decode latent representations to gene expression."""
        z_concat = np.concatenate([z_x, z_t, z_tx], axis=-1)
        x = z_concat @ self.W_dec + self.b_dec
        return x

    def get_z_tx_compositional(self, p1, p2):
        """Ground truth compositional z_tx for double KO."""
        z_tx_1 = self.z_tx_ground[p1]
        z_tx_2 = self.z_tx_ground[p2]
        if self.pathway_assignment[p1] != self.pathway_assignment[p2]:
            # Cross-pathway: additive
            return z_tx_1 + z_tx_2
        else:
            # Within-pathway: multiplicative (element-wise)
            return z_tx_1 * z_tx_2

    def generate(self, cell_type, perturbation_id, n_cells=None):
        """Generate single-cell data for a condition."""
        if n_cells is None:
            n_cells = self.n_cells

        # z_x: cell-type specific
        z_x = (np.random.randn(n_cells, self.z_dim) * self.cell_type_stds[cell_type]
               + self.cell_type_means[cell_type])

        # z_t: perturbation embedding (same for all cells in condition)
        z_t = np.tile(self.z_t_ground[perturbation_id], (n_cells, 1))

        # z_tx: interaction effect (INVARIANT across cell types)
        z_tx = np.tile(self.z_tx_ground[perturbation_id], (n_cells, 1))

        # Decode
        x = self._decoder(z_x, z_t, z_tx) + np.random.randn(n_cells, self.n_genes) * 0.3

        return x.astype(np.float32), z_x, z_t, z_tx

    def generate_double_ko(self, cell_type, p1, p2, n_cells=None):
        """Generate double-KO data."""
        if n_cells is None:
            n_cells = self.n_cells

        z_x = (np.random.randn(n_cells, self.z_dim) * self.cell_type_stds[cell_type]
               + self.cell_type_means[cell_type])

        z_t = np.tile(self.z_t_ground[p1] + self.z_t_ground[p2], (n_cells, 1))

        z_tx = np.tile(self.get_z_tx_compositional(p1, p2), (n_cells, 1))

        x = self._decoder(z_x, z_t, z_tx) + np.random.randn(n_cells, self.n_genes) * 0.3

        return x.astype(np.float32), z_x, z_t, z_tx


# ============================================================
# 2. FCR Model
# ============================================================

class FCREncoder(nn.Module):
    """Factorized Causal Representation encoder."""

    def __init__(self, n_genes, n_perturbations, z_dim, n_cell_types):
        super().__init__()
        self.z_dim = z_dim

        # Shared encoder backbone
        self.x_encoder = nn.Sequential(
            nn.Linear(n_genes, 128), nn.ReLU(),
            nn.Linear(128, 64), nn.ReLU(),
        )

        # Perturbation embedding
        self.pert_emb = nn.Embedding(n_perturbations + 10, z_dim)  # +10 for double-KO indices

        # Cell type embedding (for z_x)
        self.cell_type_emb = nn.Embedding(n_cell_types, z_dim)

        # Separate heads for z_x, z_t, z_tx
        # z_x: cell type is input (z_x SHOULD depend on cell type)
        # z_t: perturbation embedding is input
        # z_tx: NO cell type input (we want z_tx to be invariant!)
        self.z_x_head = nn.Linear(64 + n_cell_types, z_dim * 2)
        self.z_t_head = nn.Linear(64 + z_dim, z_dim * 2)
        self.z_tx_head = nn.Linear(64 + z_dim, z_dim * 2)  # NO cell type input

    def forward(self, x, pert_id, cell_type_onehot):
        h = self.x_encoder(x)
        z_t_input = self.pert_emb(pert_id)

        z_x_params = self.z_x_head(torch.cat([h, cell_type_onehot], dim=-1))
        z_x_mean, z_x_logvar = z_x_params[:, :self.z_dim], z_x_params[:, self.z_dim:]

        z_t_params = self.z_t_head(torch.cat([h, z_t_input], dim=-1))
        z_t_mean, z_t_logvar = z_t_params[:, :self.z_dim], z_t_params[:, self.z_dim:]

        z_tx_params = self.z_tx_head(torch.cat([h, z_t_input], dim=-1))
        z_tx_mean, z_tx_logvar = z_tx_params[:, :self.z_dim], z_tx_params[:, self.z_dim:]

        return (z_x_mean, z_x_logvar), (z_t_mean, z_t_logvar), (z_tx_mean, z_tx_logvar)


class FCRDecoder(nn.Module):
    """Decoder from z_x, z_t, z_tx to gene expression."""

    def __init__(self, z_dim, n_genes):
        super().__init__()
        self.decoder = nn.Sequential(
            nn.Linear(3 * z_dim, 128), nn.ReLU(),
            nn.Linear(128, n_genes),
        )

    def forward(self, z_x, z_t, z_tx):
        z = torch.cat([z_x, z_t, z_tx], dim=-1)
        return self.decoder(z)


def reparameterize(mean, logvar):
    std = torch.exp(0.5 * logvar)
    eps = torch.randn_like(std)
    return mean + eps * std


def vae_loss(x_recon, x, z_mean, z_logvar, beta=1.0):
    """ELBO loss with beta-VAE weighting."""
    recon_loss = F.mse_loss(x_recon, x, reduction='sum')
    kl = -0.5 * torch.sum(1 + z_logvar - z_mean.pow(2) - z_logvar.exp())
    return recon_loss + beta * kl


def icm_regularizer(z_tx_mean, z_tx_logvar, cell_type_onehot, cell_types):
    """
    ICM regularization: z_tx should be invariant to cell type.

    Since z_tx_head has no cell type input, we use MMD to align
    z_tx distributions across cell types.
    """
    mmd_loss = torch.tensor(0.0, device=z_tx_mean.device)
    unique_types = torch.unique(cell_types)
    if len(unique_types) > 1:
        for i in range(len(unique_types)):
            for j in range(i + 1, len(unique_types)):
                mask_i = (cell_types == unique_types[i])
                mask_j = (cell_types == unique_types[j])
                z_i = z_tx_mean[mask_i]
                z_j = z_tx_mean[mask_j]
                # Linear MMD (mean alignment)
                mmd_loss += (z_i.mean(0) - z_j.mean(0)).pow(2).sum()
                # RBF MMD
                n_sample = min(50, z_i.shape[0], z_j.shape[0])
                if n_sample > 5:
                    z_i_sub = z_i[:n_sample]
                    z_j_sub = z_j[:n_sample]
                    sigma = 1.0
                    xx = torch.exp(-torch.cdist(z_i_sub, z_i_sub).pow(2) / (2 * sigma)).mean()
                    yy = torch.exp(-torch.cdist(z_j_sub, z_j_sub).pow(2) / (2 * sigma)).mean()
                    xy = torch.exp(-torch.cdist(z_i_sub, z_j_sub).pow(2) / (2 * sigma)).mean()
                    mmd_loss += xx + yy - 2 * xy

    return mmd_loss


# ============================================================
# 3. Training
# ============================================================

def train_fcr(data, use_icm=False, n_epochs=100, lr=1e-3, beta=1.0, icm_weight=1.0):
    """Train FCR model with or without ICM regularization."""

    n_genes = data.n_genes
    n_perturbations = data.n_perturbations
    n_cell_types = data.n_cell_types
    z_dim = data.z_dim

    encoder = FCREncoder(n_genes, n_perturbations, z_dim, n_cell_types)
    decoder = FCRDecoder(z_dim, n_genes)

    optimizer = torch.optim.Adam(
        list(encoder.parameters()) + list(decoder.parameters()), lr=lr
    )

    # Generate training data: all single-KO conditions across cell types
    train_data = []
    for ct in range(n_cell_types):
        for p in range(n_perturbations):
            x, z_x, z_t, z_tx = data.generate(ct, p, n_cells=100)
            train_data.append((x, p, ct, z_tx))

    # Prepare tensors
    all_x = np.concatenate([d[0] for d in train_data])
    all_pert = np.concatenate([np.full(d[0].shape[0], d[1], dtype=np.int64) for d in train_data])
    all_ct = np.concatenate([np.full(d[0].shape[0], d[2], dtype=np.int64) for d in train_data])
    all_z_tx_gt = np.concatenate([np.tile(d[3][0], (d[0].shape[0], 1)) for d in train_data])

    x_t = torch.FloatTensor(all_x)
    pert_t = torch.LongTensor(all_pert)
    ct_t = torch.LongTensor(all_ct)
    ct_onehot = F.one_hot(ct_t, n_cell_types).float()

    z_tx_gt = torch.FloatTensor(all_z_tx_gt)

    dataset = torch.utils.data.TensorDataset(x_t, pert_t, ct_t, ct_onehot, z_tx_gt)
    loader = torch.utils.data.DataLoader(dataset, batch_size=256, shuffle=True)

    losses = []
    for epoch in range(n_epochs):
        epoch_loss = 0
        for batch_x, batch_pert, batch_ct, batch_ct_oh, batch_z_tx_gt in loader:
            optimizer.zero_grad()

            (z_x_m, z_x_lv), (z_t_m, z_t_lv), (z_tx_m, z_tx_lv) = \
                encoder(batch_x, batch_pert, batch_ct_oh)

            z_x = reparameterize(z_x_m, z_x_lv)
            z_t = reparameterize(z_t_m, z_t_lv)
            z_tx = reparameterize(z_tx_m, z_tx_lv)

            x_recon = decoder(z_x, z_t, z_tx)

            loss = vae_loss(x_recon, batch_x, z_x_m, z_x_lv, beta)
            # Add KL for z_t and z_tx (reconstruction already counted once)
            kl_t = -0.5 * torch.sum(1 + z_t_lv - z_t_m.pow(2) - z_t_lv.exp())
            kl_tx = -0.5 * torch.sum(1 + z_tx_lv - z_tx_m.pow(2) - z_tx_lv.exp())
            loss += beta * (kl_t + kl_tx)

            if use_icm:
                icm_loss = icm_regularizer(z_tx_m, z_tx_lv, batch_ct_oh, batch_ct)
                loss = loss + icm_weight * icm_loss

            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        losses.append(epoch_loss / len(loader))

        if (epoch + 1) % 20 == 0:
            print(f"  Epoch {epoch+1}/{n_epochs}: loss={losses[-1]:.2f}")

    return encoder, decoder, losses


# ============================================================
# 4. Evaluation
# ============================================================

def evaluate_invariance(encoder, data, n_cell_types, n_perturbations, z_dim):
    """RQ1: Is z_tx invariant across cell types?"""
    encoder.eval()
    results = {}

    for p in range(n_perturbations):
        z_tx_per_ct = []
        for ct in range(n_cell_types):
            x, _, _, _ = data.generate(ct, p, n_cells=200)
            x_t = torch.FloatTensor(x)
            pert_t = torch.full((200,), p, dtype=torch.long)
            ct_oh = F.one_hot(torch.full((200,), ct, dtype=torch.long), n_cell_types).float()

            with torch.no_grad():
                (_, _), (_, _), (z_tx_m, _) = encoder(x_t, pert_t, ct_oh)
            z_tx_per_ct.append(z_tx_m.mean(0).numpy())

        # Correlation between z_tx means across cell types
        if n_cell_types >= 2:
            corr = np.corrcoef(z_tx_per_ct[0], z_tx_per_ct[1])[0, 1]
            results[p] = corr

    mean_corr = np.mean(list(results.values()))
    return mean_corr, results


def evaluate_compositionality(encoder, data, n_cell_types, z_dim):
    """RQ2: Can z_tx compose to predict double-KO effects?
    Evaluation IN THE LATENT SPACE (not through decoder), since decoder
    hasn't seen composed z_tx during training.
    """
    encoder.eval()
    results = {'cross_pathway': [], 'within_pathway': []}

    for ct in range(n_cell_types):
        for p1, p2 in data.double_ko_pairs:
            # Get z_tx for each single KO
            x1, _, _, _ = data.generate(ct, p1, n_cells=200)
            x2, _, _, _ = data.generate(ct, p2, n_cells=200)

            x1_t = torch.FloatTensor(x1)
            x2_t = torch.FloatTensor(x2)
            ct_oh = F.one_hot(torch.full((200,), ct, dtype=torch.long), n_cell_types).float()

            with torch.no_grad():
                (_, _), (_, _), (z_tx1_m, _) = encoder(x1_t, torch.full((200,), p1, dtype=torch.long), ct_oh)
                (_, _), (_, _), (z_tx2_m, _) = encoder(x2_t, torch.full((200,), p2, dtype=torch.long), ct_oh)

            z_tx1 = z_tx1_m.mean(0).numpy()
            z_tx2 = z_tx2_m.mean(0).numpy()

            # Use KNOWN ground truth z_tx for double-KO (synthetic data advantage)
            z_tx_double_gt = data.get_z_tx_compositional(p1, p2)

            # Compose: try additive and multiplicative
            z_tx_add = z_tx1 + z_tx2
            z_tx_mul = z_tx1 * z_tx2

            # Evaluate: correlation between composed z_tx and encoder-inferred double-KO z_tx
            corr_add = np.corrcoef(z_tx_add, z_tx_double_gt)[0, 1]
            corr_mul = np.corrcoef(z_tx_mul, z_tx_double_gt)[0, 1]

            # Also compute cosine similarity
            def cosine_sim(a, b):
                return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-8)
            cos_add = cosine_sim(z_tx_add, z_tx_double_gt)
            cos_mul = cosine_sim(z_tx_mul, z_tx_double_gt)

            same_pathway = data.pathway_assignment[p1] == data.pathway_assignment[p2]
            key = 'within_pathway' if same_pathway else 'cross_pathway'
            results[key].append({
                'pair': (p1, p2),
                'same_pathway': same_pathway,
                'corr_additive': float(corr_add),
                'corr_multiplicative': float(corr_mul),
                'cos_additive': float(cos_add),
                'cos_multiplicative': float(cos_mul),
                'best_corr': max(float(corr_add), float(corr_mul)),
            })

    # Summary
    summary = {}
    for key in ['cross_pathway', 'within_pathway']:
        if results[key]:
            corr_add = np.mean([r['corr_additive'] for r in results[key]])
            corr_mul = np.mean([r['corr_multiplicative'] for r in results[key]])
            cos_add = np.mean([r['cos_additive'] for r in results[key]])
            cos_mul = np.mean([r['cos_multiplicative'] for r in results[key]])
            best_corr = np.mean([r['best_corr'] for r in results[key]])
            summary[key] = {
                'corr_additive': corr_add, 'corr_multiplicative': corr_mul,
                'cos_additive': cos_add, 'cos_multiplicative': cos_mul,
                'best_corr': best_corr,
            }

    return summary, results


def evaluate_transfer(encoder, data, z_dim):
    """RQ3: Zero-shot cross-cell-type transfer.
    Evaluate in latent space: does z_tx from source cell type match z_tx from target?
    """
    encoder.eval()

    source_ct, target_ct = 0, 1
    results = []

    for p in range(data.n_perturbations):
        # Source: get z_tx from cell type 0
        x_src, _, _, _ = data.generate(source_ct, p, n_cells=200)
        x_src_t = torch.FloatTensor(x_src)
        ct_oh_src = F.one_hot(torch.full((200,), source_ct, dtype=torch.long), data.n_cell_types).float()

        with torch.no_grad():
            (_, _), (_, _), (z_tx_m_src, _) = encoder(x_src_t, torch.full((200,), p, dtype=torch.long), ct_oh_src)

        z_tx_src = z_tx_m_src.mean(0).numpy()

        # Target: get z_tx from cell type 1 (ground truth for comparison)
        x_tgt, _, _, _ = data.generate(target_ct, p, n_cells=200)
        x_tgt_t = torch.FloatTensor(x_tgt)
        ct_oh_tgt = F.one_hot(torch.full((200,), target_ct, dtype=torch.long), data.n_cell_types).float()

        with torch.no_grad():
            (_, _), (_, _), (z_tx_m_tgt, _) = encoder(x_tgt_t, torch.full((200,), p, dtype=torch.long), ct_oh_tgt)

        z_tx_tgt = z_tx_m_tgt.mean(0).numpy()

        # Correlation: does source z_tx align with target z_tx?
        corr = np.corrcoef(z_tx_src, z_tx_tgt)[0, 1]

        # Cosine similarity
        cos_sim = np.dot(z_tx_src, z_tx_tgt) / (np.linalg.norm(z_tx_src) * np.linalg.norm(z_tx_tgt) + 1e-8)

        results.append({'perturbation': p, 'corr': float(corr), 'cosine': float(cos_sim)})

    mean_corr = np.mean([r['corr'] for r in results])
    mean_cos = np.mean([r['cosine'] for r in results])

    return mean_corr, mean_cos, results


# ============================================================
# 5. Main Experiment
# ============================================================

def main():
    print("=" * 70)
    print("Phase 1: Synthetic FCR+ICM Validation")
    print("=" * 70)

    # Generate synthetic data
    print("\n[1] Generating synthetic data...")
    data = SyntheticPerturbationData(
        n_genes=50, n_cell_types=2, n_perturbations=10,
        n_cells_per_condition=200, z_dim=8
    )
    print(f"  Genes: {data.n_genes}, Cell types: {data.n_cell_types}")
    print(f"  Perturbations: {data.n_perturbations} ({data.n_pathways} pathways)")
    print(f"  Double-KO pairs: {len(data.double_ko_pairs)}")
    for p1, p2 in data.double_ko_pairs:
        same = data.pathway_assignment[p1] == data.pathway_assignment[p2]
        print(f"    ({p1},{p2}): {'same' if same else 'cross'}-pathway")

    # Train FCR without ICM
    print("\n[2] Training FCR (no ICM)...")
    enc_no_icm, dec_no_icm, losses_no = train_fcr(
        data, use_icm=False, n_epochs=150, lr=1e-3, beta=0.5
    )

    # Train FCR with ICM
    print("\n[3] Training FCR + ICM...")
    enc_icm, dec_icm, losses_icm = train_fcr(
        data, use_icm=True, n_epochs=150, lr=1e-3, beta=0.5, icm_weight=10.0
    )

    # ============================================================
    # RQ1: Invariance
    # ============================================================
    print("\n" + "=" * 70)
    print("RQ1: z_tx Invariance Across Cell Types")
    print("=" * 70)

    corr_no_icm, detail_no = evaluate_invariance(enc_no_icm, data, 2, 10, 8)
    corr_icm, detail_icm = evaluate_invariance(enc_icm, data, 2, 10, 8)

    print(f"\n  FCR (no ICM): mean z_tx cross-cell correlation = {corr_no_icm:.4f}")
    print(f"  FCR + ICM:    mean z_tx cross-cell correlation = {corr_icm:.4f}")
    print(f"  Delta: {corr_icm - corr_no_icm:+.4f}")

    print(f"\n  Per-perturbation correlation:")
    print(f"  {'Pert':>6} {'No ICM':>10} {'+ICM':>10} {'Delta':>10}")
    print(f"  {'-'*36}")
    for p in range(10):
        d_no = detail_no.get(p, 0)
        d_icm = detail_icm.get(p, 0)
        print(f"  {p:>6} {d_no:>10.4f} {d_icm:>10.4f} {d_icm-d_no:>+10.4f}")

    rq1_pass = corr_icm > corr_no_icm + 0.1 and corr_icm > 0.5
    print(f"\n  RQ1 verdict: {'PASS' if rq1_pass else 'WEAK/FAIL'} "
          f"(ICM improves invariance by {corr_icm - corr_no_icm:+.4f})")

    # ============================================================
    # RQ2: Compositionality
    # ============================================================
    print("\n" + "=" * 70)
    print("RQ2: Compositional Prediction from Single-KO z_tx")
    print("=" * 70)

    summary_no, detail_no_comp = evaluate_compositionality(
        enc_no_icm, data, 2, 8)
    summary_icm, detail_icm_comp = evaluate_compositionality(
        enc_icm, data, 2, 8)

    for label, summary in [("FCR (no ICM)", summary_no), ("FCR + ICM", summary_icm)]:
        print(f"\n  {label}:")
        for key in ['cross_pathway', 'within_pathway']:
            if key in summary:
                s = summary[key]
                print(f"    {key}: corr_add={s['corr_additive']:.4f}, "
                      f"corr_mul={s['corr_multiplicative']:.4f}, "
                      f"cos_add={s['cos_additive']:.4f}, "
                      f"cos_mul={s['cos_multiplicative']:.4f}, "
                      f"best_corr={s['best_corr']:.4f}")

    rq2_pass = False
    if 'cross_pathway' in summary_icm:
        rq2_pass = summary_icm['cross_pathway']['best_corr'] > 0.5
    print(f"\n  RQ2 verdict: {'PASS' if rq2_pass else 'WEAK/FAIL'}")

    # ============================================================
    # RQ3: Zero-shot Transfer
    # ============================================================
    print("\n" + "=" * 70)
    print("RQ3: Zero-shot Cross-Cell-Type Transfer")
    print("=" * 70)

    corr_no, cos_no, detail_no_t = evaluate_transfer(
        enc_no_icm, data, 8)
    corr_icm, cos_icm, detail_icm_t = evaluate_transfer(
        enc_icm, data, 8)

    print(f"\n  FCR (no ICM): transfer corr={corr_no:.4f}, cos={cos_no:.4f}")
    print(f"  FCR + ICM:    transfer corr={corr_icm:.4f}, cos={cos_icm:.4f}")
    print(f"  Transfer improvement: {corr_icm - corr_no:+.4f}")

    rq3_pass = corr_icm > corr_no + 0.1 and corr_icm > 0.5
    print(f"\n  RQ3 verdict: {'PASS' if rq3_pass else 'WEAK/FAIL'}")

    # ============================================================
    # Overall Summary
    # ============================================================
    print("\n" + "=" * 70)
    print("PHASE 1 SUMMARY")
    print("=" * 70)

    print(f"\n  RQ1 (Invariance):       {'PASS' if rq1_pass else 'FAIL'} "
          f"(corr: {corr_no_icm:.3f} -> {corr_icm:.3f})")
    print(f"  RQ2 (Compositionality): {'PASS' if rq2_pass else 'FAIL'} "
          f"(best cross-path corr: {summary_icm.get('cross_pathway', {}).get('best_corr', 0):.3f})")
    print(f"  RQ3 (Zero-shot):        {'PASS' if rq3_pass else 'FAIL'} "
          f"(transfer corr: {corr_no:.3f} -> {corr_icm:.3f})")

    n_pass = sum([rq1_pass, rq2_pass, rq3_pass])
    print(f"\n  Overall: {n_pass}/3 RQs passed")

    if n_pass >= 2:
        print("\n  >>> SIGNAL DETECTED -- proceed to Phase 2 (real data)")
    elif n_pass == 1:
        print("\n  >>> WEAK SIGNAL -- refine ICM weight/architecture before Phase 2")
    else:
        print("\n  >>> NO SIGNAL -- hypothesis may need fundamental revision")

    return {
        "rq1_invariance": {"no_icm": corr_no_icm, "icm": corr_icm, "pass": rq1_pass},
        "rq2_composition": {"icm_best": summary_icm.get('cross_pathway', {}).get('best_corr', 0), "pass": rq2_pass},
        "rq3_transfer": {"no_icm": corr_no, "icm": corr_icm, "pass": rq3_pass},
        "n_pass": n_pass,
    }


if __name__ == "__main__":
    results = main()
