"""
Phase 2: Real Data Feasibility Check
=====================================
Test FCR+ICM on real Perturb-seq data:
  - Norman 2019: RQ2 (compositional prediction with double-KO data)
  - Replogle 2022: RQ1 (invariance) + RQ3 (cross-cell-type transfer)

Key differences from Phase 1:
  - No ground truth z_tx -> evaluate at gene expression level
  - Much higher dimensionality (5000+ genes) -> use PCA/HVG
  - Real noise, batch effects, missing conditions
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal
from sklearn.metrics import r2_score
from sklearn.decomposition import PCA
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

torch.manual_seed(42)
np.random.seed(42)

# ============================================================
# 1. Data Loading & Preprocessing
# ============================================================

def load_norman(data_path='outputs/analysis/run_04/data/gears_data'):
    """Load Norman 2019 dataset via GEARS."""
    from gears import PertData
    pd = PertData(data_path)
    pd.load(data_name='norman')
    return pd.adata


def preprocess_adata(adata, n_top_genes=500, n_pcs=50):
    """Preprocess AnnData for FCR training."""
    # Filter and normalize
    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat')
    adata = adata[:, adata.var.highly_variable].copy()

    # Store expression as float32
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray().astype(np.float32)
    else:
        X = adata.X.astype(np.float32)

    return adata, X


def get_perturbation_pairs(adata):
    """Extract single and double perturbation info."""
    conditions = adata.obs['condition'].unique()
    ctrl_key = None
    for c in conditions:
        if c == 'ctrl':
            ctrl_key = c
            break

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

    # For double perturbations, find matching single perturbations
    valid_double = []
    single_set = set(single_perts)
    for dp in double_perts:
        p1, p2 = dp.split('+')
        # Check if both single perturbations exist
        s1_exists = (p1 in single_set or f'{p1}+ctrl' in single_set or f'ctrl+{p1}' in single_set)
        s2_exists = (p2 in single_set or f'{p2}+ctrl' in single_set or f'ctrl+{p2}' in single_set)
        if s1_exists and s2_exists:
            valid_double.append(dp)

    return single_perts, valid_double, ctrl_key


# ============================================================
# 2. FCR Model (same as Phase 1, adapted for real data)
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
        self.z_tx_head = nn.Linear(128 + z_dim, z_dim * 2)  # NO cell type

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
    def __init__(self, z_dim, n_genes):
        super().__init__()
        self.decoder = nn.Sequential(
            nn.Linear(3 * z_dim, 256), nn.ReLU(),
            nn.Linear(256, n_genes),
        )

    def forward(self, z_x, z_t, z_tx):
        z = torch.cat([z_x, z_t, z_tx], dim=-1)
        return self.decoder(z)


def reparameterize(mean, logvar):
    std = torch.exp(0.5 * logvar)
    eps = torch.randn_like(std)
    return mean + eps * std


def icm_regularizer(z_tx_mean, z_tx_logvar, cell_type_onehot, cell_types):
    """MMD-based ICM regularizer."""
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

def train_fcr_real(X_train, pert_ids_train, ct_ids_train, n_perturbations,
                   n_cell_types, z_dim=16, use_icm=False, n_epochs=100,
                   lr=1e-3, beta=0.5, icm_weight=10.0, batch_size=512):
    """Train FCR on real data."""

    n_genes = X_train.shape[1]
    encoder = FCREncoder(n_genes, n_perturbations, z_dim, n_cell_types)
    decoder = FCRDecoder(z_dim, n_genes)

    optimizer = torch.optim.Adam(
        list(encoder.parameters()) + list(decoder.parameters()), lr=lr
    )

    x_t = torch.FloatTensor(X_train)
    pert_t = torch.LongTensor(pert_ids_train)
    ct_t = torch.LongTensor(ct_ids_train)
    ct_onehot = F.one_hot(ct_t, n_cell_types).float()

    dataset = torch.utils.data.TensorDataset(x_t, pert_t, ct_t, ct_onehot)
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=True)

    losses = []
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
            kl_x = -0.5 * torch.sum(1 + z_x_lv - z_x_m.pow(2) - z_x_lv.exp())
            kl_t = -0.5 * torch.sum(1 + z_t_lv - z_t_m.pow(2) - z_t_lv.exp())
            kl_tx = -0.5 * torch.sum(1 + z_tx_lv - z_tx_m.pow(2) - z_tx_lv.exp())

            loss = recon_loss + beta * (kl_x + kl_t + kl_tx)

            if use_icm and n_cell_types > 1:
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
# 4. Evaluation on Norman 2019 (RQ2: Compositionality)
# ============================================================

def evaluate_compositionality_real(encoder, decoder, adata, X, pert_id_map,
                                    n_cell_types, z_dim, n_samples=100):
    """
    RQ2 on real data: Can z_tx compose to predict double-KO gene expression?
    Evaluate at gene expression level (we have real double-KO data).
    """
    encoder.eval()
    decoder.eval()

    conditions = adata.obs['condition'].values
    unique_conditions = adata.obs['condition'].unique()

    # Identify double perturbations with matching singles
    results = []

    for dp in unique_conditions:
        if '+' not in str(dp) or 'ctrl' in str(dp):
            continue

        p1, p2 = dp.split('+')
        # Find single perturbation keys
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

        # Get real double-KO expression
        mask_double = adata.obs['condition'] == dp
        n_avail_double = mask_double.sum()
        if n_avail_double < 10:
            continue

        n_use = min(n_samples, n_avail_double)
        x_double = X[mask_double][:n_use]
        x_double_mean = x_double.mean(axis=0)

        # Get z_tx for single perturbations
        ct = 0  # Norman is single cell type
        skip = False
        for s_name, store_key in [(s1, 'z_tx_1'), (s2, 'z_tx_2')]:
            mask_s = adata.obs['condition'] == s_name
            n_avail_s = mask_s.sum()
            if n_avail_s < 10:
                skip = True
                break
            n_use_s = min(n_samples, n_avail_s)
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

        # Encode double-KO to get z_x covariate
        x_d_t = torch.FloatTensor(x_double[:n_use])
        ct_oh_d = F.one_hot(torch.full((n_use,), ct, dtype=torch.long), n_cell_types).float()
        pert_id_d = torch.full((n_use,), pert_id_map[s1], dtype=torch.long)

        with torch.no_grad():
            (z_x_m, _), (z_t_m, _), (_, _) = encoder(x_d_t, pert_id_d, ct_oh_d)

        z_x = z_x_m.mean(0).unsqueeze(0)
        z_t = z_t_m.mean(0).unsqueeze(0)

        # Compose: try additive
        z_tx_add = (z_tx_1 + z_tx_2).unsqueeze(0)
        # Compose: try multiplicative
        z_tx_mul = (z_tx_1 * z_tx_2).unsqueeze(0)

        # Decode both compositions
        with torch.no_grad():
            x_pred_add = decoder(z_x.expand(1, -1), z_t.expand(1, -1), z_tx_add)
            x_pred_mul = decoder(z_x.expand(1, -1), z_t.expand(1, -1), z_tx_mul)

        pred_add = x_pred_add[0].numpy()
        pred_mul = x_pred_mul[0].numpy()

        # Evaluate against mean double-KO expression
        r2_add = r2_score(x_double_mean, pred_add)
        r2_mul = r2_score(x_double_mean, pred_mul)

        # Also compute correlation
        corr_add = np.corrcoef(x_double_mean, pred_add)[0, 1]
        corr_mul = np.corrcoef(x_double_mean, pred_mul)[0, 1]

        results.append({
            'double_pert': dp,
            'single_1': s1, 'single_2': s2,
            'r2_additive': float(r2_add),
            'r2_multiplicative': float(r2_mul),
            'corr_additive': float(corr_add),
            'corr_multiplicative': float(corr_mul),
            'best_r2': max(float(r2_add), float(r2_mul)),
            'best_corr': max(float(corr_add), float(corr_mul)),
        })

    return results


# ============================================================
# 5. Evaluation on Replogle 2022 (RQ1 + RQ3)
# ============================================================

def evaluate_invariance_real(encoder, adata, X, pert_id_map, n_cell_types, z_dim):
    """RQ1: Is z_tx invariant across cell types?"""
    encoder.eval()

    conditions = adata.obs['condition'].unique()
    cell_types = adata.obs['cell_type'].unique()

    if len(cell_types) < 2:
        print("  WARNING: Only one cell type, cannot test invariance")
        return None

    ct0, ct1 = cell_types[0], cell_types[1]
    ct0_idx = 0
    ct1_idx = 1

    results = {}
    for cond in conditions:
        if cond == 'ctrl' or cond not in pert_id_map:
            continue

        z_tx_per_ct = []
        for ct, ct_idx in [(ct0, ct0_idx), (ct1, ct1_idx)]:
            mask = (adata.obs['condition'] == cond) & (adata.obs['cell_type'] == ct)
            if mask.sum() < 10:
                z_tx_per_ct.append(None)
                continue

            x_ct = X[mask][:100]
            x_t = torch.FloatTensor(x_ct)
            ct_oh = F.one_hot(torch.full((100,), ct_idx, dtype=torch.long), n_cell_types).float()
            pert_t = torch.full((100,), pert_id_map[cond], dtype=torch.long)

            with torch.no_grad():
                (_, _), (_, _), (z_tx_m, _) = encoder(x_t, pert_t, ct_oh)
            z_tx_per_ct.append(z_tx_m.mean(0).numpy())

        if z_tx_per_ct[0] is not None and z_tx_per_ct[1] is not None:
            corr = np.corrcoef(z_tx_per_ct[0], z_tx_per_ct[1])[0, 1]
            results[cond] = corr

    if results:
        return np.mean(list(results.values())), results
    return None, {}


# ============================================================
# 6. Main Experiment
# ============================================================

def main():
    print("=" * 70)
    print("Phase 2: Real Data Feasibility Check")
    print("=" * 70)

    # ----------------------------------------------------------
    # Part A: Norman 2019 (RQ2: Compositionality)
    # ----------------------------------------------------------
    print("\n[Part A] Norman 2019 -- RQ2: Compositional Prediction")
    print("-" * 50)

    print("  Loading Norman 2019 data...")
    adata_norman = load_norman()
    print(f"  Raw shape: {adata_norman.shape}")

    adata_norman, X_norman = preprocess_adata(adata_norman, n_top_genes=500)
    n_genes = X_norman.shape[1]
    print(f"  After HVG filter: {adata_norman.shape}")

    single_perts, double_perts, ctrl_key = get_perturbation_pairs(adata_norman)
    print(f"  Single perturbations: {len(single_perts)}")
    print(f"  Double perturbations (valid): {len(double_perts)}")

    # Build perturbation ID mapping
    all_perts = sorted(set(single_perts + double_perts + [ctrl_key]))
    pert_id_map = {p: i for i, p in enumerate(all_perts)}

    # Prepare training data (single perturbations + ctrl)
    train_conditions = [c for c in single_perts if c != ctrl_key] + [ctrl_key]
    # Also include ctrl+X as single perturbations
    train_mask = adata_norman.obs['condition'].isin(train_conditions)
    adata_train = adata_norman[train_mask].copy()
    X_train = X_norman[train_mask.values]
    pert_ids = [pert_id_map[c] for c in adata_train.obs['condition'].values]
    ct_ids = np.zeros(len(X_train), dtype=np.int64)  # Single cell type

    print(f"  Training cells: {len(X_train)}")
    print(f"  Training conditions: {len(train_conditions)}")

    # Train FCR without ICM
    print("\n  Training FCR (no ICM)...")
    enc_no, dec_no, _ = train_fcr_real(
        X_train, pert_ids, ct_ids, len(all_perts), 1,
        z_dim=16, use_icm=False, n_epochs=100, lr=1e-3, beta=0.5
    )

    # Train FCR with ICM (though single cell type, ICM should not change much)
    print("\n  Training FCR + ICM...")
    enc_icm, dec_icm, _ = train_fcr_real(
        X_train, pert_ids, ct_ids, len(all_perts), 1,
        z_dim=16, use_icm=True, n_epochs=100, lr=1e-3, beta=0.5, icm_weight=10.0
    )

    # Evaluate RQ2
    print("\n  Evaluating RQ2 (compositionality)...")
    results_no = evaluate_compositionality_real(
        enc_no, dec_no, adata_norman, X_norman, pert_id_map, 1, 16)
    results_icm = evaluate_compositionality_real(
        enc_icm, dec_icm, adata_norman, X_norman, pert_id_map, 1, 16)

    print(f"\n  FCR (no ICM): {len(results_no)} double-KO pairs tested")
    if results_no:
        mean_corr_no = np.mean([r['best_corr'] for r in results_no])
        mean_r2_no = np.mean([r['best_r2'] for r in results_no])
        print(f"    Best corr (mean): {mean_corr_no:.4f}")
        print(f"    Best R2 (mean): {mean_r2_no:.4f}")

    print(f"\n  FCR + ICM: {len(results_icm)} double-KO pairs tested")
    if results_icm:
        mean_corr_icm = np.mean([r['best_corr'] for r in results_icm])
        mean_r2_icm = np.mean([r['best_r2'] for r in results_icm])
        print(f"    Best corr (mean): {mean_corr_icm:.4f}")
        print(f"    Best R2 (mean): {mean_r2_icm:.4f}")

    # ----------------------------------------------------------
    # Part B: Replogle 2022 (RQ1 + RQ3)
    # ----------------------------------------------------------
    print("\n\n[Part B] Replogle 2022 -- RQ1 + RQ3")
    print("-" * 50)

    try:
        from gears import PertData
        pd_rep = PertData('outputs/analysis/run_04/data/gears_data')
        pd_rep.load(data_name='replogle_k562_essential')
        adata_replogle = pd_rep.adata

        print(f"  Raw shape: {adata_replogle.shape}")
        print(f"  Cell types: {adata_replogle.obs['cell_type'].unique()}")

        # If only K562, try replogle_rpe1_essential too
        if len(adata_replogle.obs['cell_type'].unique()) < 2:
            print("  Loading RPE1 data for cross-cell-type test...")
            pd_rpe1 = PertData('outputs/analysis/run_04/data/gears_data')
            pd_rpe1.load(data_name='replogle_rpe1_essential')
            adata_rpe1 = pd_rpe1.adata
            print(f"  RPE1 shape: {adata_rpe1.shape}")

            # Combine K562 and RPE1
            # Need shared perturbations
            k562_perts = set(adata_replogle.obs['condition'].unique())
            rpe1_perts = set(adata_rpe1.obs['condition'].unique())
            shared = k562_perts & rpe1_perts
            print(f"  Shared perturbations: {len(shared)}")

            if len(shared) > 5:
                # Filter to shared perturbations
                mask_k562 = adata_replogle.obs['condition'].isin(shared)
                mask_rpe1 = adata_rpe1.obs['condition'].isin(shared)

                adata_k562 = adata_replogle[mask_k562].copy()
                adata_rpe1_f = adata_rpe1[mask_rpe1].copy()

                # Set cell type
                adata_k562.obs['cell_type'] = 'K562'
                adata_rpe1_f.obs['cell_type'] = 'RPE1'

                # Concatenate (need shared genes)
                shared_genes = list(set(adata_k562.var_names) & set(adata_rpe1_f.var_names))
                print(f"  Shared genes: {len(shared_genes)}")

                if len(shared_genes) > 100:
                    adata_combined = adata_k562[:, shared_genes].concatenate(
                        adata_rpe1_f[:, shared_genes], batch_key='cell_type')
                    adata_combined, X_combined = preprocess_adata(adata_combined, n_top_genes=500)
                    print(f"  Combined shape: {adata_combined.shape}")

                    # Build pert ID map
                    all_perts_rep = sorted(adata_combined.obs['condition'].unique())
                    pert_id_map_rep = {p: i for i, p in enumerate(all_perts_rep)}

                    # Cell type IDs
                    ct_map = {'K562': 0, 'RPE1': 1}
                    ct_ids_rep = [ct_map.get(ct, 0) for ct in adata_combined.obs['cell_type'].values]
                    pert_ids_rep = [pert_id_map_rep.get(c, 0) for c in adata_combined.obs['condition'].values]

                    # Train with ICM
                    print("\n  Training FCR + ICM on combined K562+RPE1...")
                    enc_rep, dec_rep, _ = train_fcr_real(
                        X_combined, pert_ids_rep, ct_ids_rep, len(all_perts_rep), 2,
                        z_dim=16, use_icm=True, n_epochs=100, lr=1e-3, beta=0.5, icm_weight=10.0
                    )

                    # RQ1: Invariance
                    print("\n  Evaluating RQ1 (invariance)...")
                    inv_result, inv_detail = evaluate_invariance_real(
                        enc_rep, adata_combined, X_combined, pert_id_map_rep, 2, 16)
                    if inv_result is not None:
                        print(f"  Mean z_tx cross-cell correlation: {inv_result:.4f}")
                    else:
                        print("  RQ1: Could not evaluate (insufficient cross-cell data)")

                    # RQ3: Transfer
                    print("\n  Evaluating RQ3 (zero-shot transfer)...")
                    # Use latent-space evaluation like Phase 1
                    source_ct, target_ct = 0, 1
                    cell_types_rep = adata_combined.obs['cell_type'].values
                    transfer_corrs = []

                    for cond in list(shared)[:20]:
                        if cond == 'ctrl' or cond not in pert_id_map_rep:
                            continue

                        # Source (K562)
                        mask_s = (adata_combined.obs['condition'] == cond) & (cell_types_rep == 'K562')
                        # Target (RPE1)
                        mask_t = (adata_combined.obs['condition'] == cond) & (cell_types_rep == 'RPE1')

                        if mask_s.sum() < 10 or mask_t.sum() < 10:
                            continue

                        x_s = X_combined[mask_s][:100]
                        x_t = X_combined[mask_t][:100]

                        x_s_t = torch.FloatTensor(x_s)
                        x_t_t = torch.FloatTensor(x_t)

                        n_ct = 2
                        ct_oh_s = F.one_hot(torch.full((100,), 0, dtype=torch.long), n_ct).float()
                        ct_oh_t = F.one_hot(torch.full((100,), 1, dtype=torch.long), n_ct).float()
                        pert_s = torch.full((100,), pert_id_map_rep[cond], dtype=torch.long)
                        pert_t = torch.full((100,), pert_id_map_rep[cond], dtype=torch.long)

                        with torch.no_grad():
                            (_, _), (_, _), (z_tx_s, _) = enc_rep(x_s_t, pert_s, ct_oh_s)
                            (_, _), (_, _), (z_tx_t, _) = enc_rep(x_t_t, pert_t, ct_oh_t)

                        corr = np.corrcoef(z_tx_s.mean(0).numpy(), z_tx_t.mean(0).numpy())[0, 1]
                        transfer_corrs.append(corr)

                    if transfer_corrs:
                        mean_transfer_corr = np.mean(transfer_corrs)
                        print(f"  Transfer z_tx correlation (K562->RPE1): {mean_transfer_corr:.4f}")
                    else:
                        print("  RQ3: Could not evaluate (insufficient shared perturbations)")
                else:
                    print("  Not enough shared genes for cross-cell-type test")
            else:
                print("  Not enough shared perturbations for cross-cell-type test")
        else:
            # Already has multiple cell types in one dataset
            adata_replogle, X_replogle = preprocess_adata(adata_replogle, n_top_genes=500)
            print(f"  After HVG: {adata_replogle.shape}")

    except Exception as e:
        print(f"  Replogle evaluation failed: {e}")
        import traceback
        traceback.print_exc()

    # ----------------------------------------------------------
    # Summary
    # ----------------------------------------------------------
    print("\n" + "=" * 70)
    print("PHASE 2 SUMMARY")
    print("=" * 70)

    if results_icm:
        print(f"\n  Norman RQ2 (compositionality):")
        print(f"    Tested {len(results_icm)} double-KO pairs")
        print(f"    Mean best corr: {np.mean([r['best_corr'] for r in results_icm]):.4f}")
        print(f"    Mean best R2: {np.mean([r['best_r2'] for r in results_icm]):.4f}")
        # Show top 5
        top5 = sorted(results_icm, key=lambda x: x['best_corr'], reverse=True)[:5]
        print(f"    Top 5 pairs:")
        for r in top5:
            print(f"      {r['double_pert']}: corr_add={r['corr_additive']:.3f}, "
                  f"corr_mul={r['corr_multiplicative']:.3f}")


if __name__ == "__main__":
    results = main()
