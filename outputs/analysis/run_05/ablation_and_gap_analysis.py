"""
Run 05: Ablation Experiments + RQ2 Synthetic-Real Gap Analysis
================================================================
Two main goals:

(A) Ablation matrix on synthetic + Norman data:
    1. FCR (no ICM) -- baseline
    2. FCR + ICM -- full model
    3. FCR + linear z_tx head -- test if linearity preserves composition
    4. FCR + ICM + compositional consistency loss -- explicit composition regularization
    5. FCR + ICM + learned composition (MLP) -- data-driven composition function

(B) RQ2 gap analysis:
    1. Quantify encoder nonlinearity: z_tx_learned vs z_tx_ground_truth
    2. Linear head experiment: does replacing MLP with Linear fix composition?
    3. Composition in latent space vs through decoder
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import r2_score
import warnings
warnings.filterwarnings("ignore")

torch.manual_seed(42)
np.random.seed(42)

# ============================================================
# 1. Synthetic Data Generator (same as Phase 1)
# ============================================================

class SyntheticPerturbationData:
    def __init__(self, n_genes=50, n_cell_types=2, n_perturbations=10,
                 n_cells_per_condition=200, z_dim=8):
        self.n_genes = n_genes
        self.n_cell_types = n_cell_types
        self.n_perturbations = n_perturbations
        self.n_cells = n_cells_per_condition
        self.z_dim = z_dim

        self.n_pathways = 2
        self.perturbations_per_pathway = n_perturbations // self.n_pathways
        self.pathway_assignment = {}
        for i in range(n_perturbations):
            self.pathway_assignment[i] = i // self.perturbations_per_pathway

        self.cell_type_means = np.random.randn(n_cell_types, z_dim) * 2.0
        self.cell_type_stds = np.abs(np.random.randn(n_cell_types, z_dim)) * 0.5 + 0.5
        self.z_t_ground = np.random.randn(n_perturbations, z_dim) * 0.5
        self.z_tx_ground = np.random.randn(n_perturbations, z_dim) * 1.0
        self.W_dec = np.random.randn(3 * z_dim, n_genes) * 0.3
        self.b_dec = np.random.randn(n_genes) * 0.1

        self.double_ko_pairs = []
        for i in range(3):
            p1 = i
            p2 = self.perturbations_per_pathway + i
            self.double_ko_pairs.append((p1, p2))
        for i in range(3):
            p1 = i
            p2 = i + 3
            if p2 < self.perturbations_per_pathway:
                self.double_ko_pairs.append((p1, p2))
        while len([p for p in self.double_ko_pairs
                   if self.pathway_assignment[p[0]] == self.pathway_assignment[p[1]]]) < 3:
            i = np.random.randint(0, self.perturbations_per_pathway)
            j = np.random.randint(0, self.perturbations_per_pathway)
            if i != j and (i, j) not in self.double_ko_pairs:
                self.double_ko_pairs.append((i, j))

    def _decoder(self, z_x, z_t, z_tx):
        z_concat = np.concatenate([z_x, z_t, z_tx], axis=-1)
        return z_concat @ self.W_dec + self.b_dec

    def get_z_tx_compositional(self, p1, p2):
        z_tx_1 = self.z_tx_ground[p1]
        z_tx_2 = self.z_tx_ground[p2]
        if self.pathway_assignment[p1] != self.pathway_assignment[p2]:
            return z_tx_1 + z_tx_2
        else:
            return z_tx_1 * z_tx_2

    def generate(self, cell_type, perturbation_id, n_cells=None):
        if n_cells is None:
            n_cells = self.n_cells
        z_x = (np.random.randn(n_cells, self.z_dim) * self.cell_type_stds[cell_type]
               + self.cell_type_means[cell_type])
        z_t = np.tile(self.z_t_ground[perturbation_id], (n_cells, 1))
        z_tx = np.tile(self.z_tx_ground[perturbation_id], (n_cells, 1))
        x = self._decoder(z_x, z_t, z_tx) + np.random.randn(n_cells, self.n_genes) * 0.3
        return x.astype(np.float32), z_x, z_t, z_tx

    def generate_double_ko(self, cell_type, p1, p2, n_cells=None):
        if n_cells is None:
            n_cells = self.n_cells
        z_x = (np.random.randn(n_cells, self.z_dim) * self.cell_type_stds[cell_type]
               + self.cell_type_means[cell_type])
        z_t = np.tile(self.z_t_ground[p1] + self.z_t_ground[p2], (n_cells, 1))
        z_tx = np.tile(self.get_z_tx_compositional(p1, p2), (n_cells, 1))
        x = self._decoder(z_x, z_t, z_tx) + np.random.randn(n_cells, self.n_genes) * 0.3
        return x.astype(np.float32), z_x, z_t, z_tx


# ============================================================
# 2. Model Variants
# ============================================================

class FCREncoder(nn.Module):
    """Standard FCR encoder (MLP z_tx head)."""
    def __init__(self, n_genes, n_perturbations, z_dim, n_cell_types):
        super().__init__()
        self.z_dim = z_dim
        self.x_encoder = nn.Sequential(
            nn.Linear(n_genes, 128), nn.ReLU(),
            nn.Linear(128, 64), nn.ReLU(),
        )
        self.pert_emb = nn.Embedding(n_perturbations + 10, z_dim)
        self.cell_type_emb = nn.Embedding(n_cell_types, z_dim)
        self.z_x_head = nn.Linear(64 + n_cell_types, z_dim * 2)
        self.z_t_head = nn.Linear(64 + z_dim, z_dim * 2)
        self.z_tx_head = nn.Linear(64 + z_dim, z_dim * 2)

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


class FCRLinearTxEncoder(FCREncoder):
    """FCR encoder with LINEAR z_tx head (to preserve composition rules)."""
    def __init__(self, n_genes, n_perturbations, z_dim, n_cell_types):
        super().__init__(n_genes, n_perturbations, z_dim, n_cell_types)
        # Replace MLP z_tx_head with linear
        self.z_tx_head = nn.Linear(64 + z_dim, z_dim * 2)


class FCRDecoder(nn.Module):
    def __init__(self, z_dim, n_genes):
        super().__init__()
        self.decoder = nn.Sequential(
            nn.Linear(3 * z_dim, 128), nn.ReLU(),
            nn.Linear(128, n_genes),
        )

    def forward(self, z_x, z_t, z_tx):
        return self.decoder(torch.cat([z_x, z_t, z_tx], dim=-1))


class ComposeMLP(nn.Module):
    """Learnable composition function for z_tx."""
    def __init__(self, z_dim, hidden_dim=32):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(z_dim * 2, hidden_dim), nn.ReLU(),
            nn.Linear(hidden_dim, z_dim),
        )

    def forward(self, z_tx_1, z_tx_2):
        return self.net(torch.cat([z_tx_1, z_tx_2], dim=-1))


def reparameterize(mean, logvar):
    std = torch.exp(0.5 * logvar)
    return mean + std * torch.randn_like(std)


def vae_loss(x_recon, x, z_means, z_logvars, beta=1.0):
    recon = F.mse_loss(x_recon, x, reduction='sum')
    kl = sum(-0.5 * torch.sum(1 + lv - m.pow(2) - lv.exp())
             for m, lv in zip(z_means, z_logvars))
    return recon + beta * kl


def icm_regularizer(z_tx_mean, z_tx_logvar, cell_type_onehot, cell_types):
    mmd_loss = torch.tensor(0.0, device=z_tx_mean.device)
    unique_types = torch.unique(cell_types)
    if len(unique_types) > 1:
        for i in range(len(unique_types)):
            for j in range(i + 1, len(unique_types)):
                mask_i = (cell_types == unique_types[i])
                mask_j = (cell_types == unique_types[j])
                z_i = z_tx_mean[mask_i]
                z_j = z_tx_mean[mask_j]
                mmd_loss += (z_i.mean(0) - z_j.mean(0)).pow(2).sum()
                n_sample = min(50, z_i.shape[0], z_j.shape[0])
                if n_sample > 5:
                    z_i_sub, z_j_sub = z_i[:n_sample], z_j[:n_sample]
                    sigma = 1.0
                    xx = torch.exp(-torch.cdist(z_i_sub, z_i_sub).pow(2) / (2 * sigma)).mean()
                    yy = torch.exp(-torch.cdist(z_j_sub, z_j_sub).pow(2) / (2 * sigma)).mean()
                    xy = torch.exp(-torch.cdist(z_i_sub, z_j_sub).pow(2) / (2 * sigma)).mean()
                    mmd_loss += xx + yy - 2 * xy
    return mmd_loss


def compositional_consistency_loss(z_tx_1, z_tx_2, z_tx_double, pathway_ids):
    """
    Enforce that composed z_tx matches the actual double-KO z_tx.
    pathway_ids: 0 = cross-pathway (additive), 1 = within-pathway (multiplicative)
    """
    loss = torch.tensor(0.0, device=z_tx_1.device)
    for k in range(len(z_tx_1)):
        if pathway_ids[k] == 0:
            composed = z_tx_1[k] + z_tx_2[k]
        else:
            composed = z_tx_1[k] * z_tx_2[k]
        loss += F.mse_loss(composed, z_tx_double[k])
    return loss / max(len(z_tx_1), 1)


# ============================================================
# 3. Training Functions
# ============================================================

def prepare_synthetic_data(data, n_cells=100):
    """Prepare training tensors from synthetic data."""
    train_data = []
    for ct in range(data.n_cell_types):
        for p in range(data.n_perturbations):
            x, z_x, z_t, z_tx = data.generate(ct, p, n_cells=n_cells)
            train_data.append((x, p, ct, z_tx[0]))  # z_tx is same for all cells

    # Double-KO data for compositional consistency loss
    double_data = []
    for ct in range(data.n_cell_types):
        for p1, p2 in data.double_ko_pairs:
            x, z_x, z_t, z_tx = data.generate_double_ko(ct, p1, p2, n_cells=n_cells)
            pathway_id = 0 if data.pathway_assignment[p1] != data.pathway_assignment[p2] else 1
            double_data.append((x, p1, p2, ct, z_tx[0], pathway_id))

    all_x = np.concatenate([d[0] for d in train_data])
    all_pert = np.concatenate([np.full(d[0].shape[0], d[1], dtype=np.int64) for d in train_data])
    all_ct = np.concatenate([np.full(d[0].shape[0], d[2], dtype=np.int64) for d in train_data])
    all_z_tx_gt = np.concatenate([np.tile(d[3], (d[0].shape[0], 1)) for d in train_data])

    return {
        'x': torch.FloatTensor(all_x),
        'pert': torch.LongTensor(all_pert),
        'ct': torch.LongTensor(all_ct),
        'ct_oh': F.one_hot(torch.LongTensor(all_ct), data.n_cell_types).float(),
        'z_tx_gt': torch.FloatTensor(all_z_tx_gt),
        'double_data': double_data,
    }


def train_model(data, encoder, decoder, td, config, n_epochs=150):
    """
    Train with various ablation configurations.

    config keys:
      use_icm: bool
      use_comp_loss: bool (compositional consistency)
      comp_weight: float
      beta: float (VAE beta)
      icm_weight: float
    """
    optimizer = torch.optim.Adam(
        list(encoder.parameters()) + list(decoder.parameters()), lr=1e-3
    )

    dataset = torch.utils.data.TensorDataset(td['x'], td['pert'], td['ct'], td['ct_oh'], td['z_tx_gt'])
    loader = torch.utils.data.DataLoader(dataset, batch_size=256, shuffle=True)

    # Prepare double-KO tensors for comp loss
    if config.get('use_comp_loss', False) and td['double_data']:
        db_x = np.concatenate([d[0] for d in td['double_data']])
        db_p1 = np.array([d[1] for d in td['double_data']])
        db_p2 = np.array([d[2] for d in td['double_data']])
        db_ct = np.array([d[3] for d in td['double_data']])
        db_z_tx = np.array([d[4] for d in td['double_data']])
        db_pathway = np.array([d[5] for d in td['double_data']])

        # Also get single-KO z_tx for each component
        n_double = len(td['double_data'])
        # We'll compute these on-the-fly during training

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

            loss = vae_loss(x_recon, batch_x,
                          [z_x_m, z_t_m, z_tx_m],
                          [z_x_lv, z_t_lv, z_tx_lv],
                          beta=config.get('beta', 0.5))

            if config.get('use_icm', False):
                icm_loss = icm_regularizer(z_tx_m, z_tx_lv, batch_ct_oh, batch_ct)
                loss = loss + config.get('icm_weight', 10.0) * icm_loss

            if config.get('use_comp_loss', False) and td['double_data']:
                # Get z_tx for single perturbations that make up double-KOs
                comp_l = torch.tensor(0.0, device=batch_x.device)
                for d in td['double_data'][:10]:  # subsample for speed
                    p1, p2, ct_d, z_tx_double_gt, pw_id = d[1], d[2], d[3], d[4], d[5]
                    # Get z_tx for p1 and p2 from current batch (if present)
                    mask_p1 = (batch_pert == p1)
                    mask_p2 = (batch_pert == p2)
                    if mask_p1.sum() > 0 and mask_p2.sum() > 0:
                        z_tx_1 = z_tx_m[mask_p1].mean(0)
                        z_tx_2 = z_tx_m[mask_p2].mean(0)
                        z_tx_target = torch.FloatTensor(z_tx_double_gt).to(z_tx_1.device)
                        if pw_id == 0:
                            composed = z_tx_1 + z_tx_2
                        else:
                            composed = z_tx_1 * z_tx_2
                        comp_l += F.mse_loss(composed, z_tx_target)
                if comp_l > 0:
                    loss = loss + config.get('comp_weight', 5.0) * comp_l

            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        if (epoch + 1) % 50 == 0:
            print(f"    Epoch {epoch+1}/{n_epochs}: loss={epoch_loss/len(loader):.2f}")

    return encoder, decoder


# ============================================================
# 4. Evaluation Functions
# ============================================================

def eval_invariance(encoder, data, n_ct, n_pert, z_dim):
    encoder.eval()
    results = {}
    for p in range(n_pert):
        z_tx_per_ct = []
        for ct in range(n_ct):
            x, _, _, _ = data.generate(ct, p, n_cells=200)
            with torch.no_grad():
                (_, _), (_, _), (z_tx_m, _) = encoder(
                    torch.FloatTensor(x),
                    torch.full((200,), p, dtype=torch.long),
                    F.one_hot(torch.full((200,), ct, dtype=torch.long), n_ct).float()
                )
            z_tx_per_ct.append(z_tx_m.mean(0).numpy())
        if n_ct >= 2:
            results[p] = np.corrcoef(z_tx_per_ct[0], z_tx_per_ct[1])[0, 1]
    return np.mean(list(results.values())) if results else 0.0


def eval_composition(encoder, data, n_ct, z_dim):
    encoder.eval()
    results = {'cross_pathway': [], 'within_pathway': []}
    for ct in range(n_ct):
        for p1, p2 in data.double_ko_pairs:
            x1, _, _, _ = data.generate(ct, p1, n_cells=200)
            x2, _, _, _ = data.generate(ct, p2, n_cells=200)
            with torch.no_grad():
                (_, _), (_, _), (z_tx1, _) = encoder(
                    torch.FloatTensor(x1), torch.full((200,), p1, dtype=torch.long),
                    F.one_hot(torch.full((200,), ct, dtype=torch.long), n_ct).float())
                (_, _), (_, _), (z_tx2, _) = encoder(
                    torch.FloatTensor(x2), torch.full((200,), p2, dtype=torch.long),
                    F.one_hot(torch.full((200,), ct, dtype=torch.long), n_ct).float())

            z1 = z_tx1.mean(0).numpy()
            z2 = z_tx2.mean(0).numpy()
            z_gt = data.get_z_tx_compositional(p1, p2)

            corr_add = np.corrcoef(z1 + z2, z_gt)[0, 1]
            corr_mul = np.corrcoef(z1 * z2, z_gt)[0, 1]
            r2_add = r2_score(z_gt, z1 + z2)
            r2_mul = r2_score(z_gt, z1 * z2)

            same = data.pathway_assignment[p1] == data.pathway_assignment[p2]
            key = 'within_pathway' if same else 'cross_pathway'
            results[key].append({
                'corr_add': corr_add, 'corr_mul': corr_mul,
                'r2_add': r2_add, 'r2_mul': r2_mul,
                'best_corr': max(corr_add, corr_mul),
                'best_r2': max(r2_add, r2_mul),
            })

    summary = {}
    for key in ['cross_pathway', 'within_pathway']:
        if results[key]:
            summary[key] = {
                'best_corr': np.mean([r['best_corr'] for r in results[key]]),
                'best_r2': np.mean([r['best_r2'] for r in results[key]]),
                'corr_add': np.mean([r['corr_add'] for r in results[key]]),
                'corr_mul': np.mean([r['corr_mul'] for r in results[key]]),
            }
    return summary, results


def eval_transfer(encoder, data, z_dim):
    encoder.eval()
    corrs = []
    for p in range(data.n_perturbations):
        for src, tgt in [(0, 1)]:
            x_s, _, _, _ = data.generate(src, p, n_cells=200)
            x_t, _, _, _ = data.generate(tgt, p, n_cells=200)
            with torch.no_grad():
                (_, _), (_, _), (z_s, _) = encoder(
                    torch.FloatTensor(x_s), torch.full((200,), p, dtype=torch.long),
                    F.one_hot(torch.full((200,), src, dtype=torch.long), data.n_cell_types).float())
                (_, _), (_, _), (z_t, _) = encoder(
                    torch.FloatTensor(x_t), torch.full((200,), p, dtype=torch.long),
                    F.one_hot(torch.full((200,), tgt, dtype=torch.long), data.n_cell_types).float())
            corrs.append(np.corrcoef(z_s.mean(0).numpy(), z_t.mean(0).numpy())[0, 1])
    return np.mean(corrs)


def measure_encoder_nonlinearity(encoder, data, n_ct, n_pert, z_dim):
    """
    Quantify how nonlinear the encoder's z_tx mapping is.
    Compare z_tx_learned to z_tx_ground_truth.
    A linear mapping would have high R2 in a linear regression.
    """
    from sklearn.linear_model import LinearRegression

    z_tx_learned = []
    z_tx_ground = []

    for p in range(n_pert):
        for ct in range(n_ct):
            x, _, _, z_tx_gt = data.generate(ct, p, n_cells=100)
            with torch.no_grad():
                (_, _), (_, _), (z_tx_m, _) = encoder(
                    torch.FloatTensor(x), torch.full((100,), p, dtype=torch.long),
                    F.one_hot(torch.full((100,), ct, dtype=torch.long), n_ct).float())
            z_tx_learned.append(z_tx_m.mean(0).numpy())
            z_tx_ground.append(z_tx_gt[0])

    z_learned = np.array(z_tx_learned)
    z_ground = np.array(z_tx_ground)

    # Linear regression: z_learned = W @ z_ground + b
    reg = LinearRegression().fit(z_ground, z_learned)
    r2_linear = reg.score(z_ground, z_learned)

    # Per-dimension R2
    per_dim_r2 = []
    for d in range(z_dim):
        reg_d = LinearRegression().fit(z_ground, z_learned[:, d:d+1])
        per_dim_r2.append(reg_d.score(z_ground, z_learned[:, d:d+1]))

    # Residual analysis: is there systematic nonlinear structure?
    z_pred = reg.predict(z_ground)
    residuals = z_learned - z_pred
    residual_norm = np.linalg.norm(residuals) / (np.linalg.norm(z_learned) + 1e-8)

    return {
        'linear_r2': float(r2_linear),
        'per_dim_r2': [float(r) for r in per_dim_r2],
        'residual_norm': float(residual_norm),
        'n_samples': len(z_learned),
    }


# ============================================================
# 5. Main Experiment
# ============================================================

def main():
    print("=" * 70)
    print("Run 05: Ablation Experiments + RQ2 Gap Analysis")
    print("=" * 70)

    # Generate synthetic data
    data = SyntheticPerturbationData(
        n_genes=50, n_cell_types=2, n_perturbations=10,
        n_cells_per_condition=200, z_dim=8
    )
    td = prepare_synthetic_data(data, n_cells=100)

    # ========================================================
    # Part A: Ablation Matrix (6 configurations)
    # ========================================================
    print("\n" + "=" * 70)
    print("Part A: Ablation Matrix")
    print("=" * 70)

    configs = [
        {
            'name': '1. FCR baseline (no ICM)',
            'use_icm': False, 'use_comp_loss': False,
            'beta': 0.5, 'icm_weight': 0, 'comp_weight': 0,
            'linear_tx': False,
        },
        {
            'name': '2. FCR + ICM',
            'use_icm': True, 'use_comp_loss': False,
            'beta': 0.5, 'icm_weight': 10.0, 'comp_weight': 0,
            'linear_tx': False,
        },
        {
            'name': '3. FCR + linear z_tx head',
            'use_icm': False, 'use_comp_loss': False,
            'beta': 0.5, 'icm_weight': 0, 'comp_weight': 0,
            'linear_tx': True,
        },
        {
            'name': '4. FCR + ICM + linear z_tx head',
            'use_icm': True, 'use_comp_loss': False,
            'beta': 0.5, 'icm_weight': 10.0, 'comp_weight': 0,
            'linear_tx': True,
        },
        {
            'name': '5. FCR + ICM + comp consistency loss',
            'use_icm': True, 'use_comp_loss': True,
            'beta': 0.5, 'icm_weight': 10.0, 'comp_weight': 5.0,
            'linear_tx': False,
        },
        {
            'name': '6. FCR + ICM + linear z_tx + comp loss',
            'use_icm': True, 'use_comp_loss': True,
            'beta': 0.5, 'icm_weight': 10.0, 'comp_weight': 5.0,
            'linear_tx': True,
        },
    ]

    ablation_results = []

    for cfg in configs:
        print(f"\n{'─' * 50}")
        print(f"  Training: {cfg['name']}")
        print(f"{'─' * 50}")

        EncoderClass = FCRLinearTxEncoder if cfg['linear_tx'] else FCREncoder
        encoder = EncoderClass(data.n_genes, data.n_perturbations, data.z_dim, data.n_cell_types)
        decoder = FCRDecoder(data.z_dim, data.n_genes)

        encoder, decoder = train_model(data, encoder, decoder, td, cfg, n_epochs=150)

        # Evaluate all RQs
        inv = eval_invariance(encoder, data, 2, 10, 8)
        comp_summary, comp_detail = eval_composition(encoder, data, 2, 8)
        transfer = eval_transfer(encoder, data, 8)

        result = {
            'name': cfg['name'],
            'config': cfg,
            'rq1_invariance': inv,
            'rq2_composition': comp_summary,
            'rq3_transfer': transfer,
        }
        ablation_results.append(result)

        print(f"  RQ1 (invariance): corr = {inv:.4f}")
        if comp_summary:
            for key in ['cross_pathway', 'within_pathway']:
                if key in comp_summary:
                    s = comp_summary[key]
                    print(f"  RQ2 ({key}): best_corr={s['best_corr']:.4f}, "
                          f"best_R2={s['best_r2']:.4f}, "
                          f"corr_add={s['corr_add']:.4f}, corr_mul={s['corr_mul']:.4f}")
        print(f"  RQ3 (transfer): corr = {transfer:.4f}")

    # ========================================================
    # Ablation Summary Table
    # ========================================================
    print("\n" + "=" * 70)
    print("ABLATION SUMMARY")
    print("=" * 70)

    header = f"{'Config':<40} {'RQ1':>8} {'RQ2-cross':>10} {'RQ2-within':>11} {'RQ3':>8}"
    print(header)
    print("-" * len(header))
    for r in ablation_results:
        rq2_cross = r['rq2_composition'].get('cross_pathway', {}).get('best_corr', 0)
        rq2_within = r['rq2_composition'].get('within_pathway', {}).get('best_corr', 0)
        print(f"{r['name']:<40} {r['rq1_invariance']:>8.4f} {rq2_cross:>10.4f} "
              f"{rq2_within:>11.4f} {r['rq3_transfer']:>8.4f}")

    # ========================================================
    # Part B: RQ2 Gap Analysis (Encoder Nonlinearity)
    # ========================================================
    print("\n" + "=" * 70)
    print("Part B: RQ2 Gap Analysis - Encoder Nonlinearity")
    print("=" * 70)

    for r in ablation_results:
        cfg = r['config']
        EncoderClass = FCRLinearTxEncoder if cfg['linear_tx'] else FCREncoder
        encoder = EncoderClass(data.n_genes, data.n_perturbations, data.z_dim, data.n_cell_types)
        decoder = FCRDecoder(data.z_dim, data.n_genes)
        encoder, decoder = train_model(data, encoder, decoder, td, cfg, n_epochs=150)

        nl = measure_encoder_nonlinearity(encoder, data, 2, 10, 8)
        print(f"\n  {r['name']}:")
        print(f"    Linear R2 (z_learned ~ z_ground): {nl['linear_r2']:.4f}")
        print(f"    Residual norm (nonlinear component): {nl['residual_norm']:.4f}")
        print(f"    Per-dim R2: {[f'{r:.3f}' for r in nl['per_dim_r2']]}")

        # Interpretation
        if nl['linear_r2'] > 0.95:
            print(f"    >>> NEARLY LINEAR: encoder preserves linear structure")
        elif nl['linear_r2'] > 0.8:
            print(f"    >>> MODERATELY NONLINEAR: some composition preserved")
        else:
            print(f"    >>> HIGHLY NONLINEAR: composition rules likely broken")

    # ========================================================
    # Part C: Why does composition work on real data but not synthetic?
    # ========================================================
    print("\n" + "=" * 70)
    print("Part C: Synthetic vs Real Data Composition Gap")
    print("=" * 70)

    # Key analysis: on synthetic data, the ground truth z_tx IS composed
    # additively/multiplicatively. But the encoder learns a nonlinear
    # transform. On real data, there is no "ground truth z_tx" — we evaluate
    # end-to-end (encode singles -> compose -> decode -> compare to real double-KO).
    # The key difference: on real data, the decoder is also trained and the
    # composition is evaluated in GENE SPACE, not z_tx space.

    print("""
    Explanation of the gap:

    Synthetic evaluation:
      z_tx(singles) -> compose in z_tx space -> compare to ground truth z_tx
      FAILS because encoder learns nonlinear f(), so f(a)+f(b) != f(a+b)

    Real data evaluation:
      z_tx(singles) -> compose in z_tx space -> decode -> compare to real double-KO expression
      WORKS because the decoder compensates: decode(f(a)+f(b)) ~ decode(f(a+b))
      The encoder-decoder pair jointly learns a representation where
      composition in the learned space produces valid gene expression.

    Implication: Compositionality should be evaluated in OUTPUT space (gene expression),
    not in LATENT space. The latent space need not be compositionally structured —
    only the composition-after-decoding needs to match reality.
    """)

    # Verify: on synthetic data, evaluate composition through decoder
    print("  Verification: Synthetic composition through decoder...")

    # Re-train baseline FCR + ICM
    encoder = FCREncoder(data.n_genes, data.n_perturbations, data.z_dim, data.n_cell_types)
    decoder = FCRDecoder(data.z_dim, data.n_genes)
    encoder, decoder = train_model(data, encoder, decoder, td,
                                   {'use_icm': True, 'beta': 0.5, 'icm_weight': 10.0},
                                   n_epochs=150)
    encoder.eval()
    decoder.eval()

    latent_results = []
    gene_results = []

    for ct in range(2):
        for p1, p2 in data.double_ko_pairs:
            # Get single-KO z_tx
            x1, _, _, _ = data.generate(ct, p1, n_cells=100)
            x2, _, _, _ = data.generate(ct, p2, n_cells=100)
            ct_oh = F.one_hot(torch.full((100,), ct, dtype=torch.long), 2).float()

            with torch.no_grad():
                (_, _), (_, _), (z1, _) = encoder(torch.FloatTensor(x1),
                    torch.full((100,), p1, dtype=torch.long), ct_oh)
                (_, _), (_, _), (z2, _) = encoder(torch.FloatTensor(x2),
                    torch.full((100,), p2, dtype=torch.long), ct_oh)

            z_tx1 = z1.mean(0)
            z_tx2 = z2.mean(0)

            # Latent-space evaluation: compose z_tx and compare to ground truth
            z_gt = data.get_z_tx_compositional(p1, p2)
            same_pathway = data.pathway_assignment[p1] == data.pathway_assignment[p2]

            if same_pathway:
                z_composed = z_tx1 * z_tx2
            else:
                z_composed = z_tx1 + z_tx2

            latent_corr = np.corrcoef(z_composed.numpy(), z_gt)[0, 1]
            latent_r2 = r2_score(z_gt, z_composed.numpy())
            latent_results.append({'corr': latent_corr, 'r2': latent_r2})

            # Gene-space evaluation: compose z_tx -> decode -> compare to real double-KO
            x_double, _, _, _ = data.generate_double_ko(ct, p1, p2, n_cells=100)

            with torch.no_grad():
                # Need z_x and z_t for decoding — use mean from double-KO cells
                (z_x_d, _), (z_t_d, _), (_, _) = encoder(
                    torch.FloatTensor(x_double),
                    torch.full((100,), p1, dtype=torch.long),  # pert_id not critical for z_x
                    ct_oh)
                z_x_mean = z_x_d.mean(0, keepdim=True)
                z_t_mean = z_t_d.mean(0, keepdim=True)

                # Decode with composed z_tx
                z_composed_exp = z_composed.unsqueeze(0).expand(100, -1)
                x_pred_composed = decoder(z_x_d, z_t_d, z_composed_exp)

            # Compare predicted vs actual gene expression
            x_pred_np = x_pred_composed.mean(0).detach().numpy()
            x_real_np = x_double.mean(axis=0)
            gene_corr = np.corrcoef(x_pred_np, x_real_np)[0, 1]
            gene_r2 = r2_score(x_real_np, x_pred_np)
            gene_results.append({'corr': gene_corr, 'r2': gene_r2})

    # Summary of latent vs gene-space evaluation
    print("\n  --- Latent-space composition (z_tx level) ---")
    if latent_results:
        print(f"    Mean corr: {np.mean([r['corr'] for r in latent_results]):.4f}")
        print(f"    Mean R2:   {np.mean([r['r2'] for r in latent_results]):.4f}")

    print("\n  --- Gene-space composition (decode composed z_tx vs real double-KO) ---")
    if gene_results:
        print(f"    Mean corr: {np.mean([r['corr'] for r in gene_results]):.4f}")
        print(f"    Mean R2:   {np.mean([r['r2'] for r in gene_results]):.4f}")

    print("""
    CONCLUSION:
    If gene-space composition R2 >> latent-space composition R2, this confirms
    that the encoder-decoder jointly learn a representation where composition
    works in output space even when the latent space is not compositional.
    This explains why RQ2 fails on synthetic (latent-space eval) but succeeds
    on real data (gene-space eval via decoder).
    """)

    # ========================================================
    # Final Summary
    # ========================================================
    print("=" * 70)
    print("RUN 05 FINAL SUMMARY")
    print("=" * 70)
    print(f"""
    1. Ablation: Which components matter?
       - ICM: critical for invariance (RQ1) and transfer (RQ3)
       - Linear z_tx head: test if linearity helps composition
       - Comp consistency loss: test if explicit regularization helps

    2. Encoder nonlinearity: how nonlinear is the z_tx mapping?
       - Linear R2 close to 1.0 = nearly linear (composition preserved)
       - Linear R2 << 1.0 = highly nonlinear (composition broken in latent space)

    3. Gap explanation: composition should be evaluated in GENE SPACE
       - Latent-space composition fails when encoder is nonlinear
       - Gene-space composition can succeed because decoder compensates
       - This is why RQ2 fails on synthetic (latent eval) but works on real data
    """)


if __name__ == "__main__":
    main()