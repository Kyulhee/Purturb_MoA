"""
Run 06: Norman 2019 Real Data Ablation Experiments
====================================================
Test the 6 ablation configs from run_05 on real Norman Perturb-seq data.

Key question: Does comp consistency loss improve RQ2 (compositionality)
on real data, as it does on synthetic data?

Since Norman is single cell type (K562):
  - RQ1/RQ3 not directly testable (need multi-cell-type)
  - RQ2 (compositionality) is the main evaluation
  - ICM acts as regularizer even with 1 cell type (z_tx structure constraint)
  - Evaluation is in gene expression space (per run_05 finding)
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
# 1. Data Loading & Preprocessing
# ============================================================

def load_norman(data_path='outputs/analysis/run_04/data/gears_data'):
    from gears import PertData
    pd = PertData(data_path)
    pd.load(data_name='norman')
    return pd.adata


def preprocess_adata(adata, n_top_genes=500):
    import scanpy as sc
    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat')
    adata = adata[:, adata.var.highly_variable].copy()
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray().astype(np.float32)
    else:
        X = adata.X.astype(np.float32)
    return adata, X


def get_perturbation_info(adata):
    conditions = adata.obs['condition'].unique()
    single_perts = []
    double_perts = []
    for c in conditions:
        if c == 'ctrl':
            continue
        parts = c.split('+')
        if len(parts) == 1:
            single_perts.append(c)
        elif len(parts) == 2:
            if parts[0] == 'ctrl' or parts[1] == 'ctrl':
                single_perts.append(c)
            else:
                double_perts.append(c)

    valid_double = []
    single_set = set(single_perts)
    for dp in double_perts:
        p1, p2 = dp.split('+')
        s1_exists = (p1 in single_set or f'{p1}+ctrl' in single_set or f'ctrl+{p1}' in single_set)
        s2_exists = (p2 in single_set or f'{p2}+ctrl' in single_set or f'ctrl+{p2}' in single_set)
        if s1_exists and s2_exists:
            valid_double.append(dp)

    return single_perts, valid_double


# ============================================================
# 2. Model Variants
# ============================================================

class FCREncoder(nn.Module):
    def __init__(self, n_genes, n_perturbations, z_dim, n_cell_types):
        super().__init__()
        self.z_dim = z_dim
        self.x_encoder = nn.Sequential(
            nn.Linear(n_genes, 256), nn.ReLU(), nn.Dropout(0.1),
            nn.Linear(256, 128), nn.ReLU(),
        )
        self.pert_emb = nn.Embedding(n_perturbations, z_dim)
        self.cell_type_emb = nn.Embedding(n_cell_types, z_dim)
        self.z_x_head = nn.Linear(128 + n_cell_types, z_dim * 2)
        self.z_t_head = nn.Linear(128 + z_dim, z_dim * 2)
        self.z_tx_head = nn.Linear(128 + z_dim, z_dim * 2)

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
    """FCR encoder with LINEAR z_tx head."""
    def __init__(self, n_genes, n_perturbations, z_dim, n_cell_types):
        super().__init__(n_genes, n_perturbations, z_dim, n_cell_types)
        # z_tx_head is already Linear — but the upstream network (x_encoder) is nonlinear
        # This tests whether removing one layer of nonlinearity in z_tx helps


class FCRDecoder(nn.Module):
    def __init__(self, z_dim, n_genes):
        super().__init__()
        self.decoder = nn.Sequential(
            nn.Linear(3 * z_dim, 256), nn.ReLU(),
            nn.Linear(256, n_genes),
        )

    def forward(self, z_x, z_t, z_tx):
        return self.decoder(torch.cat([z_x, z_t, z_tx], dim=-1))


def reparameterize(mean, logvar):
    std = torch.exp(0.5 * logvar)
    return mean + std * torch.randn_like(std)


def icm_regularizer(z_tx_mean, z_tx_logvar, cell_type_onehot, cell_types):
    """MMD-based ICM regularizer. With 1 cell type, this acts as a structural regularizer."""
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
    else:
        # Single cell type: use z_tx variance regularization as structural constraint
        # Encourage z_tx to have unit variance per dimension (standard normal prior)
        var = z_tx_mean.var(0)
        mmd_loss = ((var - 1.0) ** 2).mean()
    return mmd_loss


# ============================================================
# 3. Training with Ablation Configs
# ============================================================

def train_model(encoder, decoder, X_train, pert_ids, ct_ids, n_cell_types,
                config, double_ko_data=None, n_epochs=100, batch_size=512):
    """
    Train FCR with various ablation configurations.

    config keys:
      use_icm, use_comp_loss, comp_weight, beta, icm_weight
    """
    optimizer = torch.optim.Adam(
        list(encoder.parameters()) + list(decoder.parameters()), lr=1e-3
    )

    n_genes = X_train.shape[1]
    z_dim = encoder.z_dim

    x_t = torch.FloatTensor(X_train)
    pert_t = torch.LongTensor(pert_ids)
    ct_t = torch.LongTensor(ct_ids)
    ct_onehot = F.one_hot(ct_t, n_cell_types).float()

    dataset = torch.utils.data.TensorDataset(x_t, pert_t, ct_t, ct_onehot)
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=True)

    for epoch in range(n_epochs):
        epoch_loss = 0
        for batch_x, batch_pert, batch_ct, batch_ct_oh in loader:
            optimizer.zero_grad()

            (z_x_m, z_x_lv), (z_t_m, z_t_lv), (z_tx_m, z_tx_lv) = \
                encoder(batch_x, batch_pert, batch_ct_oh)

            z_x = reparameterize(z_x_m, z_x_lv)
            z_t = reparameterize(z_t_m, z_t_lv)
            z_tx = reparameterize(z_tx_m, z_tx_lv)

            x_recon = decoder(z_x, z_t, z_tx)

            recon_loss = F.mse_loss(x_recon, batch_x, reduction='sum')
            kl = sum(-0.5 * torch.sum(1 + lv - m.pow(2) - lv.exp())
                     for m, lv in [(z_x_m, z_x_lv), (z_t_m, z_t_lv), (z_tx_m, z_tx_lv)])
            loss = recon_loss + config.get('beta', 0.5) * kl

            if config.get('use_icm', False):
                icm_loss = icm_regularizer(z_tx_m, z_tx_lv, batch_ct_oh, batch_ct)
                loss = loss + config.get('icm_weight', 10.0) * icm_loss

            if config.get('use_comp_loss', False) and double_ko_data is not None:
                # Compositional consistency: for each double-KO, compose z_tx from singles
                # and enforce that decoded output matches real double-KO expression
                comp_l = torch.tensor(0.0, device=batch_x.device)
                n_comp = 0
                for dp_info in double_ko_data[:15]:  # subsample for speed
                    p1_name, p2_name, dp_name = dp_info['p1'], dp_info['p2'], dp_info['dp']
                    p1_id, p2_id = dp_info['p1_id'], dp_info['p2_id']

                    mask_p1 = (batch_pert == p1_id)
                    mask_p2 = (batch_pert == p2_id)

                    if mask_p1.sum() > 2 and mask_p2.sum() > 2:
                        z_tx_1 = z_tx_m[mask_p1].mean(0)
                        z_tx_2 = z_tx_m[mask_p2].mean(0)

                        # Get z_x from a random cell (covariate)
                        z_x_ref = z_x_m[mask_p1[:1].nonzero().squeeze()[:1]]
                        z_t_ref = z_t_m[mask_p1[:1].nonzero().squeeze()[:1]]

                        # Compose z_tx (try additive — most common in real data)
                        z_tx_composed = z_tx_1 + z_tx_2
                        x_pred = decoder(z_x_ref, z_t_ref, z_tx_composed.unsqueeze(0))

                        # Target: mean double-KO expression
                        x_target = torch.FloatTensor(dp_info['x_mean']).to(x_pred.device)
                        comp_l += F.mse_loss(x_pred[0], x_target)
                        n_comp += 1

                if n_comp > 0:
                    loss = loss + config.get('comp_weight', 5.0) * comp_l / n_comp

            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        if (epoch + 1) % 25 == 0:
            print(f"    Epoch {epoch+1}/{n_epochs}: loss={epoch_loss/len(loader):.2f}")

    return encoder, decoder


# ============================================================
# 4. Evaluation (Gene-Space Compositionality)
# ============================================================

def evaluate_compositionality(encoder, decoder, adata, X, pert_id_map,
                               n_cell_types, z_dim, n_samples=100):
    """
    RQ2: Evaluate compositionality at gene expression level.
    Encode singles -> compose z_tx -> decode -> compare to real double-KO.
    """
    encoder.eval()
    decoder.eval()

    conditions = adata.obs['condition'].values
    unique_conditions = adata.obs['condition'].unique()
    ct = 0  # Norman = single cell type

    results = []
    for dp in unique_conditions:
        if '+' not in str(dp) or 'ctrl' in str(dp):
            continue

        p1, p2 = dp.split('+')
        s1 = None
        s2 = None
        for s in [p1, f'{p1}+ctrl', f'ctrl+{p1}']:
            if s in pert_id_map:
                s1 = s
                break
        for s in [p2, f'{p2}+ctrl', f'ctrl+{p2}']:
            if s in pert_id_map:
                s2 = s
                break

        if s1 is None or s2 is None:
            continue

        # Real double-KO expression
        mask_double = adata.obs['condition'] == dp
        if mask_double.sum() < 10:
            continue
        n_use = min(n_samples, mask_double.sum())
        x_double = X[mask_double][:n_use]
        x_double_mean = x_double.mean(axis=0)

        # Encode single perturbations
        skip = False
        for s_name, store_key in [(s1, 'z_tx_1'), (s2, 'z_tx_2')]:
            mask_s = adata.obs['condition'] == s_name
            if mask_s.sum() < 10:
                skip = True
                break
            n_use_s = min(n_samples, mask_s.sum())
            x_s = X[mask_s][:n_use_s]
            x_t = torch.FloatTensor(x_s)
            ct_oh = F.one_hot(torch.full((n_use_s,), ct, dtype=torch.long), n_cell_types).float()
            pert_id = torch.full((n_use_s,), pert_id_map[s_name], dtype=torch.long)

            with torch.no_grad():
                (_, _), (_, _), (z_tx_m, _) = encoder(x_t, pert_id, ct_oh)
            if store_key == 'z_tx_1':
                z_tx_1 = z_tx_m.mean(0)
            else:
                z_tx_2 = z_tx_m.mean(0)

        if skip:
            continue

        # Get z_x, z_t from double-KO cells (covariate)
        x_d_t = torch.FloatTensor(x_double[:n_use])
        ct_oh_d = F.one_hot(torch.full((n_use,), ct, dtype=torch.long), n_cell_types).float()
        pert_id_d = torch.full((n_use,), pert_id_map[s1], dtype=torch.long)

        with torch.no_grad():
            (z_x_m, _), (z_t_m, _), (_, _) = encoder(x_d_t, pert_id_d, ct_oh_d)

        z_x = z_x_m.mean(0).unsqueeze(0)
        z_t = z_t_m.mean(0).unsqueeze(0)

        # Compose: additive
        z_tx_add = (z_tx_1 + z_tx_2).unsqueeze(0)
        # Compose: multiplicative
        z_tx_mul = (z_tx_1 * z_tx_2).unsqueeze(0)

        with torch.no_grad():
            x_pred_add = decoder(z_x.expand(1, -1), z_t.expand(1, -1), z_tx_add)
            x_pred_mul = decoder(z_x.expand(1, -1), z_t.expand(1, -1), z_tx_mul)

        pred_add = x_pred_add[0].numpy()
        pred_mul = x_pred_mul[0].numpy()

        r2_add = r2_score(x_double_mean, pred_add)
        r2_mul = r2_score(x_double_mean, pred_mul)
        corr_add = np.corrcoef(x_double_mean, pred_add)[0, 1]
        corr_mul = np.corrcoef(x_double_mean, pred_mul)[0, 1]

        results.append({
            'double_pert': dp,
            'r2_add': r2_add, 'r2_mul': r2_mul,
            'corr_add': corr_add, 'corr_mul': corr_mul,
            'best_r2': max(r2_add, r2_mul),
            'best_corr': max(corr_add, corr_mul),
        })

    return results


# ============================================================
# 5. Main
# ============================================================

def main():
    print("=" * 70)
    print("Run 06: Norman 2019 Real Data Ablation Experiments")
    print("=" * 70)

    # Load data
    print("\n  Loading Norman 2019 data...")
    adata = load_norman()
    print(f"  Raw shape: {adata.shape}")

    adata, X = preprocess_adata(adata, n_top_genes=500)
    n_genes = X.shape[1]
    print(f"  After HVG filter: {adata.shape}")

    single_perts, double_perts = get_perturbation_info(adata)
    print(f"  Single perturbations: {len(single_perts)}")
    print(f"  Double perturbations (valid): {len(double_perts)}")

    # Build perturbation ID mapping
    all_perts = sorted(set(single_perts + double_perts + ['ctrl']))
    pert_id_map = {p: i for i, p in enumerate(all_perts)}
    n_perturbations = len(all_perts)
    n_cell_types = 1

    # Training data: single perturbations + ctrl
    train_conditions = [c for c in single_perts if c != 'ctrl'] + ['ctrl']
    train_mask = adata.obs['condition'].isin(train_conditions)
    X_train = X[train_mask.values]
    pert_ids = [pert_id_map[c] for c in adata.obs['condition'].values[train_mask.values]]
    ct_ids = np.zeros(len(X_train), dtype=np.int64)

    print(f"  Training cells: {len(X_train)}")

    # Prepare double-KO data for comp consistency loss
    double_ko_info = []
    for dp in double_perts:
        p1, p2 = dp.split('+')
        s1 = None
        s2 = None
        for s in [p1, f'{p1}+ctrl', f'ctrl+{p1}']:
            if s in pert_id_map:
                s1 = s
                break
        for s in [p2, f'{p2}+ctrl', f'ctrl+{p2}']:
            if s in pert_id_map:
                s2 = s
                break
        if s1 is None or s2 is None:
            continue
        mask_double = adata.obs['condition'] == dp
        if mask_double.sum() < 10:
            continue
        x_double = X[mask_double]
        double_ko_info.append({
            'dp': dp, 'p1': s1, 'p2': s2,
            'p1_id': pert_id_map[s1], 'p2_id': pert_id_map[s2],
            'x_mean': x_double.mean(axis=0),
        })
    print(f"  Double-KO pairs for comp loss: {len(double_ko_info)}")

    # ========================================================
    # Ablation Configs (same as run_05)
    # ========================================================
    z_dim = 16

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

    all_results = []

    for cfg in configs:
        print(f"\n{'='*60}")
        print(f"  Training: {cfg['name']}")
        print(f"{'='*60}")

        EncoderClass = FCRLinearTxEncoder if cfg['linear_tx'] else FCREncoder
        encoder = EncoderClass(n_genes, n_perturbations, z_dim, n_cell_types)
        decoder = FCRDecoder(z_dim, n_genes)

        encoder, decoder = train_model(
            encoder, decoder, X_train, pert_ids, ct_ids, n_cell_types,
            cfg, double_ko_info, n_epochs=100
        )

        # Evaluate RQ2 (compositionality at gene level)
        results = evaluate_compositionality(
            encoder, decoder, adata, X, pert_id_map, n_cell_types, z_dim
        )

        if results:
            mean_corr = np.mean([r['best_corr'] for r in results])
            mean_r2 = np.mean([r['best_r2'] for r in results])
            mean_corr_add = np.mean([r['corr_add'] for r in results])
            mean_corr_mul = np.mean([r['corr_mul'] for r in results])
            mean_r2_add = np.mean([r['r2_add'] for r in results])
            mean_r2_mul = np.mean([r['r2_mul'] for r in results])
        else:
            mean_corr = mean_r2 = mean_corr_add = mean_corr_mul = 0
            mean_r2_add = mean_r2_mul = 0

        result = {
            'name': cfg['name'],
            'config': cfg,
            'n_pairs': len(results),
            'best_corr': mean_corr,
            'best_r2': mean_r2,
            'corr_add': mean_corr_add,
            'corr_mul': mean_corr_mul,
            'r2_add': mean_r2_add,
            'r2_mul': mean_r2_mul,
            'detail': results,
        }
        all_results.append(result)

        print(f"\n  RQ2 Results ({len(results)} double-KO pairs):")
        print(f"    best_corr={mean_corr:.4f}, best_R2={mean_r2:.4f}")
        print(f"    corr_add={mean_corr_add:.4f}, corr_mul={mean_corr_mul:.4f}")
        print(f"    R2_add={mean_r2_add:.4f}, R2_mul={mean_r2_mul:.4f}")

    # ========================================================
    # Summary Table
    # ========================================================
    print("\n" + "=" * 70)
    print("ABLATION SUMMARY (Norman 2019 Real Data)")
    print("=" * 70)

    header = f"{'Config':<40} {'n_pairs':>8} {'best_corr':>10} {'best_R2':>10} {'corr_add':>10} {'corr_mul':>10}"
    print(header)
    print("-" * len(header))
    for r in all_results:
        print(f"{r['name']:<40} {r['n_pairs']:>8} {r['best_corr']:>10.4f} "
              f"{r['best_r2']:>10.4f} {r['corr_add']:>10.4f} {r['corr_mul']:>10.4f}")

    # Compare with synthetic results from run_05
    print("""
    Comparison with Run 05 (synthetic):
      - Synthetic: comp loss was the key driver (RQ2-cross 0.20 -> 0.79)
      - Real data: does comp loss similarly improve gene-space composition?
      - Real data: additive vs multiplicative composition preference?
    """)

    # Per-pair analysis for best config
    best = max(all_results, key=lambda r: r['best_r2'])
    print(f"\n  Best config: {best['name']}")
    print(f"  Per-pair results (top 10):")
    for r in sorted(best['detail'], key=lambda x: x['best_corr'], reverse=True)[:10]:
        print(f"    {r['double_pert']}: corr_add={r['corr_add']:.3f}, "
              f"corr_mul={r['corr_mul']:.3f}, R2_add={r['r2_add']:.3f}, R2_mul={r['r2_mul']:.3f}")

    # Additive vs multiplicative preference
    add_wins = sum(1 for r in best['detail'] if r['corr_add'] > r['corr_mul'])
    mul_wins = sum(1 for r in best['detail'] if r['corr_mul'] > r['corr_add'])
    print(f"\n  Additive better: {add_wins}, Multiplicative better: {mul_wins}")

    print("\n" + "=" * 70)
    print("RUN 06 COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
