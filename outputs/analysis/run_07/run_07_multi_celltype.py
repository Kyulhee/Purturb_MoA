"""
Run 07: Multi-Cell-Type FCR-ICM Experiment
============================================
RQ1 + RQ3 validation on real data using Replogle 2022 K562+RPE1.

Key fix from run_04: AnnData.concatenate(batch_key='cell_type') overwrites
cell_type values with batch indices. Fix: use batch_key='batch' instead.

Optimizations vs run_04 script:
  - Subsample per perturbation (max 200 cells) to keep memory manageable
  - Flush prints for real-time progress
  - Load h5ad directly (skip PyG) for speed

Configs:
  1. FCR baseline (no ICM)
  2. FCR + ICM
"""

import sys
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import r2_score
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

torch.manual_seed(42)
np.random.seed(42)

def flush_print(*args, **kwargs):
    print(*args, **kwargs, flush=True)


# ============================================================
# 1. Data Loading
# ============================================================

def load_and_combine(data_path='outputs/analysis/run_04/data/gears_data',
                     max_cells_per_pert=200):
    """Load K562 and RPE1, combine with shared perturbations/genes."""
    from gears import PertData

    flush_print("  Loading K562...")
    pd_k562 = PertData(data_path)
    pd_k562.load(data_name='replogle_k562_essential')
    adata_k562 = pd_k562.adata
    flush_print(f"    K562 shape: {adata_k562.shape}")

    flush_print("  Loading RPE1...")
    pd_rpe1 = PertData(data_path)
    pd_rpe1.load(data_name='replogle_rpe1_essential')
    adata_rpe1 = pd_rpe1.adata
    flush_print(f"    RPE1 shape: {adata_rpe1.shape}")

    # Shared perturbations
    k562_perts = set(adata_k562.obs['condition'].unique())
    rpe1_perts = set(adata_rpe1.obs['condition'].unique())
    shared_perts = k562_perts & rpe1_perts
    shared_perts.discard('ctrl')
    flush_print(f"  Shared perturbations (excl ctrl): {len(shared_perts)}")

    # Shared genes — intersect BEFORE any heavy processing
    shared_genes = sorted(set(adata_k562.var_names) & set(adata_rpe1.var_names))
    flush_print(f"  Shared genes: {len(shared_genes)}")

    # Filter to shared perturbations + shared genes
    keep_perts = shared_perts | {'ctrl'}
    mask_k = adata_k562.obs['condition'].isin(keep_perts)
    mask_r = adata_rpe1.obs['condition'].isin(keep_perts)

    adata_k = adata_k562[mask_k, shared_genes].copy()
    adata_r = adata_rpe1[mask_r, shared_genes].copy()
    flush_print(f"  After filtering: K562={adata_k.shape}, RPE1={adata_r.shape}")

    # Subsample per perturbation to keep memory manageable
    flush_print(f"  Subsampling (max {max_cells_per_pert} cells/pert)...")
    adata_k = subsample_per_condition(adata_k, max_cells_per_pert)
    adata_r = subsample_per_condition(adata_r, max_cells_per_pert)
    flush_print(f"  After subsampling: K562={adata_k.shape}, RPE1={adata_r.shape}")

    # Set cell type BEFORE concatenate
    adata_k.obs['cell_type'] = 'K562'
    adata_r.obs['cell_type'] = 'RPE1'

    # CRITICAL FIX: batch_key='batch', NOT 'cell_type'
    adata_combined = adata_k.concatenate(adata_r, batch_key='batch')
    flush_print(f"  Combined shape: {adata_combined.shape}")
    flush_print(f"  Cell types after concat: {adata_combined.obs['cell_type'].unique()}")

    return adata_combined, shared_perts


def subsample_per_condition(adata, max_cells):
    """Subsample to max_cells per condition."""
    indices = []
    for cond in adata.obs['condition'].unique():
        mask = adata.obs['condition'] == cond
        idx = np.where(mask)[0]
        if len(idx) > max_cells:
            idx = np.random.choice(idx, max_cells, replace=False)
        indices.extend(idx)
    return adata[sorted(indices)].copy()


def preprocess_adata(adata, n_top_genes=500):
    """Preprocess combined AnnData."""
    flush_print("  Normalizing and selecting HVGs...")
    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat')
    adata = adata[:, adata.var.highly_variable].copy()
    flush_print(f"  After HVG: {adata.shape}")

    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray().astype(np.float32)
    else:
        X = adata.X.astype(np.float32)

    return adata, X


# ============================================================
# 2. FCR Model
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
        return self.decoder(torch.cat([z_x, z_t, z_tx], dim=-1))


def reparameterize(mean, logvar):
    std = torch.exp(0.5 * logvar)
    return mean + torch.randn_like(std) * std


def icm_regularizer(z_tx_mean, cell_types):
    """MMD-based ICM regularizer — aligns z_tx across cell types."""
    mmd_loss = torch.tensor(0.0, device=z_tx_mean.device)
    unique_types = torch.unique(cell_types)
    if len(unique_types) < 2:
        return mmd_loss

    for i in range(len(unique_types)):
        for j in range(i + 1, len(unique_types)):
            mask_i = (cell_types == unique_types[i])
            mask_j = (cell_types == unique_types[j])
            z_i = z_tx_mean[mask_i]
            z_j = z_tx_mean[mask_j]

            mmd_loss += (z_i.mean(0) - z_j.mean(0)).pow(2).sum()

            n_sample = min(100, z_i.shape[0], z_j.shape[0])
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

def train_fcr(X_train, pert_ids_train, ct_ids_train, n_perturbations,
              n_cell_types, z_dim=8, use_icm=False, n_epochs=150,
              lr=1e-3, beta=0.5, icm_weight=10.0, batch_size=512):
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
            kl = -0.5 * torch.sum(1 + z_x_lv - z_x_m.pow(2) - z_x_lv.exp()) + \
                 -0.5 * torch.sum(1 + z_t_lv - z_t_m.pow(2) - z_t_lv.exp()) + \
                 -0.5 * torch.sum(1 + z_tx_lv - z_tx_m.pow(2) - z_tx_lv.exp())

            loss = recon_loss + beta * kl

            if use_icm and n_cell_types > 1:
                icm_loss = icm_regularizer(z_tx_m, batch_ct)
                loss = loss + icm_weight * icm_loss

            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        if (epoch + 1) % 30 == 0:
            flush_print(f"    Epoch {epoch+1}/{n_epochs}: loss={epoch_loss/len(loader):.1f}")

    return encoder, decoder


# ============================================================
# 4. RQ1: z_tx Invariance
# ============================================================

def evaluate_rq1(encoder, adata, X, pert_id_map, n_cell_types, z_dim, n_samples=100):
    """RQ1: Is z_tx invariant across cell types?"""
    encoder.eval()
    ct_names = sorted(adata.obs['cell_type'].unique())
    ct_map = {name: i for i, name in enumerate(ct_names)}

    results = []
    for cond in adata.obs['condition'].unique():
        if cond == 'ctrl' or cond not in pert_id_map:
            continue

        z_tx_per_ct = []
        for ct in ct_names:
            mask = (adata.obs['condition'] == cond) & (adata.obs['cell_type'] == ct)
            if mask.sum() < 10:
                z_tx_per_ct.append(None)
                continue
            n_use = min(n_samples, mask.sum())
            x_ct = X[mask][:n_use]
            x_t = torch.FloatTensor(x_ct)
            ct_oh = F.one_hot(torch.full((n_use,), ct_map[ct], dtype=torch.long), n_cell_types).float()
            pert_t = torch.full((n_use,), pert_id_map[cond], dtype=torch.long)

            with torch.no_grad():
                (_, _), (_, _), (z_tx_m, _) = encoder(x_t, pert_t, ct_oh)
            z_tx_per_ct.append(z_tx_m.mean(0).numpy())

        if z_tx_per_ct[0] is not None and z_tx_per_ct[1] is not None:
            corr = np.corrcoef(z_tx_per_ct[0], z_tx_per_ct[1])[0, 1]
            cos = np.dot(z_tx_per_ct[0], z_tx_per_ct[1]) / (
                np.linalg.norm(z_tx_per_ct[0]) * np.linalg.norm(z_tx_per_ct[1]) + 1e-8)
            results.append({'pert': cond, 'corr': float(corr), 'cosine': float(cos)})

    if results:
        return {
            'mean_corr': np.mean([r['corr'] for r in results]),
            'mean_cosine': np.mean([r['cosine'] for r in results]),
            'n_tested': len(results),
            'details': results,
        }
    return None


# ============================================================
# 5. RQ3: Zero-Shot Transfer
# ============================================================

def evaluate_rq3(encoder, decoder, adata, X, pert_id_map, n_cell_types, z_dim,
                 source_ct='K562', target_ct='RPE1', n_samples=100):
    """RQ3: Use source z_tx + target z_x to predict target expression."""
    encoder.eval()
    decoder.eval()

    ct_map = {source_ct: 0, target_ct: 1}

    results = []
    for cond in adata.obs['condition'].unique():
        if cond == 'ctrl' or cond not in pert_id_map:
            continue

        mask_s = (adata.obs['condition'] == cond) & (adata.obs['cell_type'] == source_ct)
        mask_t = (adata.obs['condition'] == cond) & (adata.obs['cell_type'] == target_ct)
        if mask_s.sum() < 10 or mask_t.sum() < 10:
            continue

        # Source z_tx
        n_s = min(n_samples, mask_s.sum())
        x_s = X[mask_s][:n_s]
        x_s_t = torch.FloatTensor(x_s)
        ct_oh_s = F.one_hot(torch.full((n_s,), 0, dtype=torch.long), n_cell_types).float()
        pert_s = torch.full((n_s,), pert_id_map[cond], dtype=torch.long)
        with torch.no_grad():
            (z_x_s_m, _), (z_t_s_m, _), (z_tx_s_m, _) = encoder(x_s_t, pert_s, ct_oh_s)

        # Target z_x, z_t
        n_t = min(n_samples, mask_t.sum())
        x_t_data = X[mask_t][:n_t]
        x_t_t = torch.FloatTensor(x_t_data)
        ct_oh_t = F.one_hot(torch.full((n_t,), 1, dtype=torch.long), n_cell_types).float()
        pert_t = torch.full((n_t,), pert_id_map[cond], dtype=torch.long)
        with torch.no_grad():
            (z_x_t_m, _), (z_t_t_m, _), (z_tx_t_m, _) = encoder(x_t_t, pert_t, ct_oh_t)

        # Transfer: source z_tx + target z_x, z_t
        with torch.no_grad():
            x_pred = decoder(
                z_x_t_m.mean(0).unsqueeze(0),
                z_t_t_m.mean(0).unsqueeze(0),
                z_tx_s_m.mean(0).unsqueeze(0),
            )

        # Oracle: target's own z_tx
        with torch.no_grad():
            x_oracle = decoder(
                z_x_t_m.mean(0).unsqueeze(0),
                z_t_t_m.mean(0).unsqueeze(0),
                z_tx_t_m.mean(0).unsqueeze(0),
            )

        actual = x_t_data.mean(axis=0)
        pred = x_pred[0].numpy()
        oracle = x_oracle[0].numpy()

        r2_t = r2_score(actual, pred)
        corr_t = np.corrcoef(actual, pred)[0, 1]
        r2_o = r2_score(actual, oracle)
        corr_o = np.corrcoef(actual, oracle)[0, 1]
        z_tx_corr = np.corrcoef(z_tx_s_m.mean(0).numpy(), z_tx_t_m.mean(0).numpy())[0, 1]

        results.append({
            'pert': cond,
            'r2_transfer': float(r2_t), 'corr_transfer': float(corr_t),
            'r2_oracle': float(r2_o), 'corr_oracle': float(corr_o),
            'z_tx_cross_corr': float(z_tx_corr),
        })

    if results:
        return {
            'mean_r2_transfer': np.mean([r['r2_transfer'] for r in results]),
            'mean_corr_transfer': np.mean([r['corr_transfer'] for r in results]),
            'mean_r2_oracle': np.mean([r['r2_oracle'] for r in results]),
            'mean_corr_oracle': np.mean([r['corr_oracle'] for r in results]),
            'mean_z_tx_cross_corr': np.mean([r['z_tx_cross_corr'] for r in results]),
            'n_tested': len(results),
            'details': results,
        }
    return None


# ============================================================
# 6. Main
# ============================================================

def main():
    flush_print("=" * 70)
    flush_print("Run 07: Multi-Cell-Type FCR-ICM (RQ1 + RQ3 on Real Data)")
    flush_print("=" * 70)

    # Load
    flush_print("\n[1] Loading Replogle K562 + RPE1...")
    adata, shared_perts = load_and_combine(max_cells_per_pert=200)
    flush_print(f"  Cell type counts: {adata.obs['cell_type'].value_counts().to_dict()}")

    # Preprocess
    flush_print("\n[2] Preprocessing...")
    adata, X = preprocess_adata(adata, n_top_genes=500)
    n_genes = X.shape[1]

    # Mappings
    ct_names = sorted(adata.obs['cell_type'].unique())
    ct_map = {name: i for i, name in enumerate(ct_names)}
    n_cell_types = len(ct_names)

    all_perts = sorted(adata.obs['condition'].unique())
    pert_id_map = {p: i for i, p in enumerate(all_perts)}
    n_perts = len(all_perts)

    ct_ids = np.array([ct_map[ct] for ct in adata.obs['cell_type'].values], dtype=np.int64)
    pert_ids = np.array([pert_id_map[c] for c in adata.obs['condition'].values], dtype=np.int64)

    for ct_name, ct_id in ct_map.items():
        flush_print(f"  {ct_name} (id={ct_id}): {(ct_ids == ct_id).sum()} cells")
    flush_print(f"  Perturbations: {n_perts}")

    # Config 1: FCR baseline
    flush_print("\n" + "=" * 50)
    flush_print("[3] Config 1: FCR baseline (no ICM)")
    enc1, dec1 = train_fcr(
        X, pert_ids, ct_ids, n_perts, n_cell_types,
        z_dim=8, use_icm=False, n_epochs=150
    )

    # Config 2: FCR + ICM
    flush_print("\n" + "=" * 50)
    flush_print("[4] Config 2: FCR + ICM")
    enc2, dec2 = train_fcr(
        X, pert_ids, ct_ids, n_perts, n_cell_types,
        z_dim=8, use_icm=True, n_epochs=150, icm_weight=10.0
    )

    # RQ1
    flush_print("\n" + "=" * 50)
    flush_print("[5] RQ1: z_tx Invariance Across Cell Types")
    rq1_no = evaluate_rq1(enc1, adata, X, pert_id_map, n_cell_types, 8)
    rq1_icm = evaluate_rq1(enc2, adata, X, pert_id_map, n_cell_types, 8)

    if rq1_no:
        flush_print(f"  FCR baseline: corr={rq1_no['mean_corr']:.4f}, cosine={rq1_no['mean_cosine']:.4f} ({rq1_no['n_tested']} perts)")
    if rq1_icm:
        flush_print(f"  FCR + ICM:    corr={rq1_icm['mean_corr']:.4f}, cosine={rq1_icm['mean_cosine']:.4f} ({rq1_icm['n_tested']} perts)")

    # RQ3
    flush_print("\n" + "=" * 50)
    flush_print("[6] RQ3: Zero-Shot Transfer (K562 -> RPE1)")
    rq3_no = evaluate_rq3(enc1, dec1, adata, X, pert_id_map, n_cell_types, 8)
    rq3_icm = evaluate_rq3(enc2, dec2, adata, X, pert_id_map, n_cell_types, 8)

    if rq3_no:
        flush_print(f"  FCR baseline: R2={rq3_no['mean_r2_transfer']:.4f}, corr={rq3_no['mean_corr_transfer']:.4f}, z_tx_corr={rq3_no['mean_z_tx_cross_corr']:.4f} ({rq3_no['n_tested']} perts)")
    if rq3_icm:
        flush_print(f"  FCR + ICM:    R2={rq3_icm['mean_r2_transfer']:.4f}, corr={rq3_icm['mean_corr_transfer']:.4f}, z_tx_corr={rq3_icm['mean_z_tx_cross_corr']:.4f} ({rq3_icm['n_tested']} perts)")

    # Summary
    flush_print("\n" + "=" * 70)
    flush_print("RUN 07 SUMMARY")
    flush_print("=" * 70)

    summary = {
        'dataset': 'Replogle 2022 (K562 + RPE1)',
        'n_shared_perts': len(shared_perts),
        'n_genes_hvg': n_genes,
        'n_cells_total': X.shape[0],
    }

    if rq1_no and rq1_icm:
        summary['rq1_no_icm_corr'] = rq1_no['mean_corr']
        summary['rq1_icm_corr'] = rq1_icm['mean_corr']
        summary['rq1_no_icm_cosine'] = rq1_no['mean_cosine']
        summary['rq1_icm_cosine'] = rq1_icm['mean_cosine']
        flush_print(f"\n  RQ1 (Invariance):")
        flush_print(f"    FCR baseline: corr={rq1_no['mean_corr']:.4f}, cosine={rq1_no['mean_cosine']:.4f}")
        flush_print(f"    FCR + ICM:    corr={rq1_icm['mean_corr']:.4f}, cosine={rq1_icm['mean_cosine']:.4f}")

    if rq3_no and rq3_icm:
        summary['rq3_no_icm_r2'] = rq3_no['mean_r2_transfer']
        summary['rq3_icm_r2'] = rq3_icm['mean_r2_transfer']
        summary['rq3_no_icm_corr'] = rq3_no['mean_corr_transfer']
        summary['rq3_icm_corr'] = rq3_icm['mean_corr_transfer']
        summary['rq3_no_icm_ztx_corr'] = rq3_no['mean_z_tx_cross_corr']
        summary['rq3_icm_ztx_corr'] = rq3_icm['mean_z_tx_cross_corr']
        flush_print(f"\n  RQ3 (Transfer K562->RPE1):")
        flush_print(f"    FCR baseline: R2={rq3_no['mean_r2_transfer']:.4f}, corr={rq3_no['mean_corr_transfer']:.4f}")
        flush_print(f"    FCR + ICM:    R2={rq3_icm['mean_r2_transfer']:.4f}, corr={rq3_icm['mean_corr_transfer']:.4f}")
        flush_print(f"    z_tx cross-corr: baseline={rq3_no['mean_z_tx_cross_corr']:.4f}, ICM={rq3_icm['mean_z_tx_cross_corr']:.4f}")

    # Save
    import json
    save_path = 'outputs/analysis/run_07/run_07_results.json'
    with open(save_path, 'w') as f:
        json.dump(summary, f, default=lambda o: float(o) if isinstance(o, np.floating) else int(o) if isinstance(o, np.integer) else o, indent=2)
    flush_print(f"\n  Results saved to {save_path}")

    # Top improvements
    if rq1_no and rq1_icm:
        no_d = {d['pert']: d['corr'] for d in rq1_no['details']}
        icm_d = {d['pert']: d['corr'] for d in rq1_icm['details']}
        improvements = [(p, no_d[p], icm_d[p], icm_d[p]-no_d[p]) for p in no_d if p in icm_d]
        improvements.sort(key=lambda x: x[3], reverse=True)
        flush_print(f"\n  Top 10 z_tx corr improvements:")
        for p, b, a, d in improvements[:10]:
            flush_print(f"    {p}: {b:.3f} -> {a:.3f} (delta={d:+.3f})")

    return summary


if __name__ == "__main__":
    main()
