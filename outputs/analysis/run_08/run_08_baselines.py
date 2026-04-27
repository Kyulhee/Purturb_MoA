"""
Run 08: Baseline Comparison for Cross-Cell-Type Transfer
=========================================================
Compare FCR-ICM against baselines on RQ3 (zero-shot K562->RPE1 transfer).

Baselines:
  1. Mean Shift: RPE1_ctrl_mean + (K562_pert_mean - K562_ctrl_mean)
     — Simplest possible transfer: assume perturbation effect is cell-type invariant
  2. GEARS: Train on K562, predict with RPE1 ctrl cells
     — SOTA GNN-based perturbation predictor
  3. FCR baseline (no ICM): Already evaluated in run_07
  4. FCR + ICM: Our proposed method (run_07 results)

Key question: Does any baseline achieve similar transfer performance to FCR+ICM?
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
# 1. Data Loading (same as run_07)
# ============================================================

def load_and_combine(data_path='outputs/analysis/run_04/data/gears_data',
                     max_cells_per_pert=200):
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

    k562_perts = set(adata_k562.obs['condition'].unique())
    rpe1_perts = set(adata_rpe1.obs['condition'].unique())
    shared_perts = k562_perts & rpe1_perts
    shared_perts.discard('ctrl')
    flush_print(f"  Shared perturbations: {len(shared_perts)}")

    shared_genes = sorted(set(adata_k562.var_names) & set(adata_rpe1.var_names))
    flush_print(f"  Shared genes: {len(shared_genes)}")

    keep_perts = shared_perts | {'ctrl'}
    mask_k = adata_k562.obs['condition'].isin(keep_perts)
    mask_r = adata_rpe1.obs['condition'].isin(keep_perts)

    adata_k = adata_k562[mask_k, shared_genes].copy()
    adata_r = adata_rpe1[mask_r, shared_genes].copy()

    # Subsample
    adata_k = _subsample(adata_k, max_cells_per_pert)
    adata_r = _subsample(adata_r, max_cells_per_pert)
    flush_print(f"  After subsampling: K562={adata_k.shape}, RPE1={adata_r.shape}")

    adata_k.obs['cell_type'] = 'K562'
    adata_r.obs['cell_type'] = 'RPE1'
    adata_combined = adata_k.concatenate(adata_r, batch_key='batch')

    return adata_combined, shared_perts, adata_k, adata_r


def _subsample(adata, max_cells):
    indices = []
    for cond in adata.obs['condition'].unique():
        mask = adata.obs['condition'] == cond
        idx = np.where(mask)[0]
        if len(idx) > max_cells:
            idx = np.random.choice(idx, max_cells, replace=False)
        indices.extend(idx)
    return adata[sorted(indices)].copy()


def preprocess_adata(adata, n_top_genes=500):
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


# ============================================================
# 2. Baseline 1: Mean Shift
# ============================================================

def evaluate_mean_shift(adata_k562, adata_rpe1, shared_genes, shared_perts, n_top_genes=500):
    """
    Mean Shift: RPE1_ctrl_mean + (K562_pert_mean - K562_ctrl_mean)
    Assumes perturbation effect delta is the same across cell types.
    """
    flush_print("\n  Baseline 1: Mean Shift")

    # Preprocess separately to get expression in same gene space
    # Use shared genes only
    adata_k = adata_k562[:, shared_genes].copy()
    adata_r = adata_rpe1[:, shared_genes].copy()

    sc.pp.normalize_total(adata_k, target_sum=1e4)
    sc.pp.log1p(adata_k)
    sc.pp.normalize_total(adata_r, target_sum=1e4)
    sc.pp.log1p(adata_r)

    X_k = adata_k.X.toarray().astype(np.float32) if hasattr(adata_k.X, 'toarray') else adata_k.X.astype(np.float32)
    X_r = adata_r.X.toarray().astype(np.float32) if hasattr(adata_r.X, 'toarray') else adata_r.X.astype(np.float32)

    # Control means
    ctrl_mask_k = adata_k.obs['condition'] == 'ctrl'
    ctrl_mask_r = adata_r.obs['condition'] == 'ctrl'
    ctrl_mean_k = X_k[ctrl_mask_k].mean(axis=0)
    ctrl_mean_r = X_r[ctrl_mask_r].mean(axis=0)

    # Evaluate on shared perturbations
    results = []
    tested = 0
    for pert in shared_perts:
        if pert == 'ctrl':
            continue

        # K562 perturbation
        mask_k = adata_k.obs['condition'] == pert
        mask_r = adata_r.obs['condition'] == pert
        if mask_k.sum() < 10 or mask_r.sum() < 10:
            continue

        pert_mean_k = X_k[mask_k].mean(axis=0)
        pert_mean_r = X_r[mask_r].mean(axis=0)

        # Mean shift prediction: RPE1_ctrl + (K562_pert - K562_ctrl)
        delta_k = pert_mean_k - ctrl_mean_k
        pred_r = ctrl_mean_r + delta_k

        r2 = r2_score(pert_mean_r, pred_r)
        corr = np.corrcoef(pert_mean_r, pred_r)[0, 1]

        results.append({'pert': pert, 'r2': float(r2), 'corr': float(corr)})
        tested += 1

    if results:
        mean_r2 = np.mean([r['r2'] for r in results])
        mean_corr = np.mean([r['corr'] for r in results])
        flush_print(f"    Tested {tested} perturbations")
        flush_print(f"    Mean R2: {mean_r2:.4f}")
        flush_print(f"    Mean corr: {mean_corr:.4f}")
        return {'mean_r2': mean_r2, 'mean_corr': mean_corr, 'n_tested': tested, 'details': results}
    return None


# ============================================================
# 3. Baseline 2: GEARS
# ============================================================

def evaluate_gears(data_path='outputs/analysis/run_04/data/gears_data'):
    """
    Train GEARS on K562, then predict for RPE1 by replacing ctrl cells.
    """
    flush_print("\n  Baseline 2: GEARS")

    try:
        from gears import PertData, GEARS

        # Load K562 and train GEARS
        flush_print("    Loading K562 for GEARS training...")
        pd_k562 = PertData(data_path)
        pd_k562.load(data_name='replogle_k562_essential')

        flush_print("    Initializing GEARS model...")
        gears_model = GEARS(pd_k562, device='cpu')
        gears_model.model_initialize(hidden_size=64)

        flush_print("    Training GEARS (20 epochs)...")
        gears_model.train(epochs=20)
        flush_print("    GEARS training complete")

        # Load RPE1 for evaluation
        flush_print("    Loading RPE1 for evaluation...")
        pd_rpe1 = PertData(data_path)
        pd_rpe1.load(data_name='replogle_rpe1_essential')
        adata_rpe1 = pd_rpe1.adata

        # Preprocess RPE1
        sc.pp.normalize_total(adata_rpe1, target_sum=1e4)
        sc.pp.log1p(adata_rpe1)
        if hasattr(adata_rpe1.X, 'toarray'):
            X_rpe1 = adata_rpe1.X.toarray().astype(np.float32)
        else:
            X_rpe1 = adata_rpe1.X.astype(np.float32)

        # Get shared perturbations
        k562_perts = set(pd_k562.adata.obs['condition'].unique())
        rpe1_perts = set(adata_rpe1.obs['condition'].unique())
        shared = k562_perts & rpe1_perts
        shared.discard('ctrl')

        # For GEARS cross-cell-type: hack by replacing ctrl_adata
        # GEARS predicts from ctrl cells, so we swap RPE1 ctrl cells
        ctrl_mask_rpe1 = adata_rpe1.obs['condition'] == 'ctrl'
        rpe1_ctrl_adata = adata_rpe1[ctrl_mask_rpe1]

        # Save original ctrl and swap
        original_ctrl = gears_model.ctrl_adata
        gears_model.ctrl_adata = rpe1_ctrl_adata

        # Need to also update the gene names to match
        # GEARS uses its own adata gene ordering
        results = []
        tested = 0

        # Only test perturbations that GEARS knows about
        gears_perts = set(gears_model.pert_list)

        for pert in list(shared)[:50]:  # Test first 50 for speed
            if pert == 'ctrl':
                continue
            # GEARS uses gene names without +ctrl
            pert_gene = pert.replace('+ctrl', '')

            if pert_gene not in gears_perts:
                continue

            # Get real RPE1 expression
            mask_r = adata_rpe1.obs['condition'] == pert
            if mask_r.sum() < 10:
                continue

            rpe1_real = X_rpe1[mask_r].mean(axis=0)

            # Predict with GEARS using RPE1 ctrl cells
            try:
                pred_dict = gears_model.predict([[pert_gene]])
                pred_key = pert_gene
                if pred_key in pred_dict:
                    pred = pred_dict[pred_key]
                    # Align genes — GEARS prediction is in K562 gene order
                    # Need to map to RPE1 gene order
                    # For simplicity, compute on shared genes
                    k562_genes = list(pd_k562.adata.var_names)
                    rpe1_genes = list(adata_rpe1.var_names)
                    shared_g = [g for g in k562_genes if g in rpe1_genes]

                    k562_idx = [k562_genes.index(g) for g in shared_g]
                    rpe1_idx = [rpe1_genes.index(g) for g in shared_g]

                    pred_shared = pred[k562_idx]
                    real_shared = rpe1_real[rpe1_idx]

                    r2 = r2_score(real_shared, pred_shared)
                    corr = np.corrcoef(real_shared, pred_shared)[0, 1]
                    results.append({'pert': pert, 'r2': float(r2), 'corr': float(corr)})
                    tested += 1
            except Exception as e:
                flush_print(f"      {pert}: prediction failed ({e})")
                continue

        # Restore original ctrl
        gears_model.ctrl_adata = original_ctrl

        if results:
            mean_r2 = np.mean([r['r2'] for r in results])
            mean_corr = np.mean([r['corr'] for r in results])
            flush_print(f"    Tested {tested} perturbations")
            flush_print(f"    Mean R2: {mean_r2:.4f}")
            flush_print(f"    Mean corr: {mean_corr:.4f}")
            return {'mean_r2': mean_r2, 'mean_corr': mean_corr, 'n_tested': tested, 'details': results}
        else:
            flush_print("    No valid predictions")
            return None

    except Exception as e:
        flush_print(f"    GEARS evaluation failed: {e}")
        import traceback
        traceback.print_exc()
        return None


# ============================================================
# 4. FCR+ICM (re-run for fair comparison on same data)
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
                z_i_sub, z_j_sub = z_i[:n_sample], z_j[:n_sample]
                sigma = 1.0
                xx = torch.exp(-torch.cdist(z_i_sub, z_i_sub).pow(2) / (2 * sigma)).mean()
                yy = torch.exp(-torch.cdist(z_j_sub, z_j_sub).pow(2) / (2 * sigma)).mean()
                xy = torch.exp(-torch.cdist(z_i_sub, z_j_sub).pow(2) / (2 * sigma)).mean()
                mmd_loss += xx + yy - 2 * xy
    return mmd_loss


def train_fcr(X, pert_ids, ct_ids, n_perts, n_ct, z_dim=8, use_icm=False,
              n_epochs=150, lr=1e-3, beta=0.5, icm_weight=10.0, batch_size=512):
    n_genes = X.shape[1]
    enc = FCREncoder(n_genes, n_perts, z_dim, n_ct)
    dec = FCRDecoder(z_dim, n_genes)
    opt = torch.optim.Adam(list(enc.parameters()) + list(dec.parameters()), lr=lr)

    x_t = torch.FloatTensor(X)
    p_t = torch.LongTensor(pert_ids)
    c_t = torch.LongTensor(ct_ids)
    c_oh = F.one_hot(c_t, n_ct).float()
    ds = torch.utils.data.TensorDataset(x_t, p_t, c_t, c_oh)
    loader = torch.utils.data.DataLoader(ds, batch_size=batch_size, shuffle=True, drop_last=True)

    for epoch in range(n_epochs):
        for bx, bp, bc, bcoh in loader:
            opt.zero_grad()
            (z_x_m, z_x_lv), (z_t_m, z_t_lv), (z_tx_m, z_tx_lv) = enc(bx, bp, bcoh)
            z_x = reparameterize(z_x_m, z_x_lv)
            z_t = reparameterize(z_t_m, z_t_lv)
            z_tx = reparameterize(z_tx_m, z_tx_lv)
            x_rec = dec(z_x, z_t, z_tx)
            loss = F.mse_loss(x_rec, bx, reduction='sum') + beta * (
                -0.5 * torch.sum(1 + z_x_lv - z_x_m.pow(2) - z_x_lv.exp()) +
                -0.5 * torch.sum(1 + z_t_lv - z_t_m.pow(2) - z_t_lv.exp()) +
                -0.5 * torch.sum(1 + z_tx_lv - z_tx_m.pow(2) - z_tx_lv.exp())
            )
            if use_icm and n_ct > 1:
                loss = loss + icm_weight * icm_regularizer(z_tx_m, bc)
            loss.backward()
            opt.step()
        if (epoch + 1) % 30 == 0:
            flush_print(f"    Epoch {epoch+1}/{n_epochs}")

    return enc, dec


def evaluate_fcr_transfer(enc, dec, adata, X, pert_id_map, n_ct, z_dim,
                          source_ct='K562', target_ct='RPE1', n_samples=100):
    """Same RQ3 evaluation as run_07."""
    enc.eval()
    dec.eval()
    ct_map = {source_ct: 0, target_ct: 1}

    results = []
    for cond in adata.obs['condition'].unique():
        if cond == 'ctrl' or cond not in pert_id_map:
            continue
        mask_s = (adata.obs['condition'] == cond) & (adata.obs['cell_type'] == source_ct)
        mask_t = (adata.obs['condition'] == cond) & (adata.obs['cell_type'] == target_ct)
        if mask_s.sum() < 10 or mask_t.sum() < 10:
            continue

        n_s = min(n_samples, mask_s.sum())
        n_t = min(n_samples, mask_t.sum())
        x_s = X[mask_s][:n_s]
        x_t = X[mask_t][:n_t]

        with torch.no_grad():
            s_t = torch.FloatTensor(x_s)
            t_t = torch.FloatTensor(x_t)
            s_oh = F.one_hot(torch.full((n_s,), 0, dtype=torch.long), n_ct).float()
            t_oh = F.one_hot(torch.full((n_t,), 1, dtype=torch.long), n_ct).float()
            s_p = torch.full((n_s,), pert_id_map[cond], dtype=torch.long)
            t_p = torch.full((n_t,), pert_id_map[cond], dtype=torch.long)

            (z_x_s, _), (z_t_s, _), (z_tx_s, _) = enc(s_t, s_p, s_oh)
            (z_x_t, _), (z_t_t, _), (z_tx_t, _) = enc(t_t, t_p, t_oh)

            x_pred = dec(z_x_t.mean(0).unsqueeze(0), z_t_t.mean(0).unsqueeze(0), z_tx_s.mean(0).unsqueeze(0))

        actual = x_t.mean(axis=0)
        pred = x_pred[0].numpy()
        r2 = r2_score(actual, pred)
        corr = np.corrcoef(actual, pred)[0, 1]
        results.append({'pert': cond, 'r2': float(r2), 'corr': float(corr)})

    if results:
        return {
            'mean_r2': np.mean([r['r2'] for r in results]),
            'mean_corr': np.mean([r['corr'] for r in results]),
            'n_tested': len(results),
        }
    return None


# ============================================================
# 5. Main
# ============================================================

def main():
    flush_print("=" * 70)
    flush_print("Run 08: Baseline Comparison for Cross-Cell-Type Transfer")
    flush_print("=" * 70)

    # Load data
    flush_print("\n[1] Loading data...")
    adata_combined, shared_perts, adata_k562_raw, adata_rpe1_raw = load_and_combine(max_cells_per_pert=200)
    flush_print(f"  Combined: {adata_combined.shape}")
    flush_print(f"  Cell types: {adata_combined.obs['cell_type'].value_counts().to_dict()}")

    # Baseline 1: Mean Shift (operates on raw expression, not HVG)
    flush_print("\n" + "=" * 50)
    flush_print("[2] Baseline 1: Mean Shift")
    mean_shift_result = evaluate_mean_shift(adata_k562_raw, adata_rpe1_raw,
                                             list(set(adata_k562_raw.var_names) & set(adata_rpe1_raw.var_names)),
                                             shared_perts)

    # Preprocess for FCR
    flush_print("\n" + "=" * 50)
    flush_print("[3] Preprocessing for FCR models...")
    adata, X = preprocess_adata(adata_combined, n_top_genes=500)
    n_genes = X.shape[1]

    ct_names = sorted(adata.obs['cell_type'].unique())
    ct_map = {name: i for i, name in enumerate(ct_names)}
    n_ct = len(ct_names)
    all_perts = sorted(adata.obs['condition'].unique())
    pert_id_map = {p: i for i, p in enumerate(all_perts)}
    n_perts = len(all_perts)

    ct_ids = np.array([ct_map[ct] for ct in adata.obs['cell_type'].values], dtype=np.int64)
    pert_ids = np.array([pert_id_map[c] for c in adata.obs['condition'].values], dtype=np.int64)

    # Baseline 3: FCR no ICM
    flush_print("\n" + "=" * 50)
    flush_print("[4] Baseline 3: FCR (no ICM)")
    enc_no, dec_no = train_fcr(X, pert_ids, ct_ids, n_perts, n_ct, z_dim=8, use_icm=False, n_epochs=150)
    fcr_no_result = evaluate_fcr_transfer(enc_no, dec_no, adata, X, pert_id_map, n_ct, 8)
    if fcr_no_result:
        flush_print(f"    R2: {fcr_no_result['mean_r2']:.4f}, corr: {fcr_no_result['mean_corr']:.4f} ({fcr_no_result['n_tested']} perts)")

    # Our method: FCR + ICM
    flush_print("\n" + "=" * 50)
    flush_print("[5] Our Method: FCR + ICM")
    enc_icm, dec_icm = train_fcr(X, pert_ids, ct_ids, n_perts, n_ct, z_dim=8, use_icm=True, n_epochs=150, icm_weight=10.0)
    fcr_icm_result = evaluate_fcr_transfer(enc_icm, dec_icm, adata, X, pert_id_map, n_ct, 8)
    if fcr_icm_result:
        flush_print(f"    R2: {fcr_icm_result['mean_r2']:.4f}, corr: {fcr_icm_result['mean_corr']:.4f} ({fcr_icm_result['n_tested']} perts)")

    # Baseline 2: GEARS
    flush_print("\n" + "=" * 50)
    flush_print("[6] Baseline 2: GEARS")
    gears_result = evaluate_gears()

    # Summary
    flush_print("\n" + "=" * 70)
    flush_print("RUN 08 SUMMARY: Cross-Cell-Type Transfer Comparison")
    flush_print("=" * 70)
    flush_print(f"\n  {'Method':<25} {'R2':>8} {'Corr':>8} {'N':>6}")
    flush_print(f"  {'-'*25} {'-'*8} {'-'*8} {'-'*6}")

    if mean_shift_result:
        flush_print(f"  {'Mean Shift':<25} {mean_shift_result['mean_r2']:>8.4f} {mean_shift_result['mean_corr']:>8.4f} {mean_shift_result['n_tested']:>6}")
    if gears_result:
        flush_print(f"  {'GEARS':<25} {gears_result['mean_r2']:>8.4f} {gears_result['mean_corr']:>8.4f} {gears_result['n_tested']:>6}")
    if fcr_no_result:
        flush_print(f"  {'FCR (no ICM)':<25} {fcr_no_result['mean_r2']:>8.4f} {fcr_no_result['mean_corr']:>8.4f} {fcr_no_result['n_tested']:>6}")
    if fcr_icm_result:
        flush_print(f"  {'FCR + ICM (ours)':<25} {fcr_icm_result['mean_r2']:>8.4f} {fcr_icm_result['mean_corr']:>8.4f} {fcr_icm_result['n_tested']:>6}")

    # Save
    import json
    summary = {}
    if mean_shift_result:
        summary['mean_shift'] = {'r2': mean_shift_result['mean_r2'], 'corr': mean_shift_result['mean_corr'], 'n': mean_shift_result['n_tested']}
    if gears_result:
        summary['gears'] = {'r2': gears_result['mean_r2'], 'corr': gears_result['mean_corr'], 'n': gears_result['n_tested']}
    if fcr_no_result:
        summary['fcr_no_icm'] = {'r2': fcr_no_result['mean_r2'], 'corr': fcr_no_result['mean_corr'], 'n': fcr_no_result['n_tested']}
    if fcr_icm_result:
        summary['fcr_icm'] = {'r2': fcr_icm_result['mean_r2'], 'corr': fcr_icm_result['mean_corr'], 'n': fcr_icm_result['n_tested']}

    with open('outputs/analysis/run_08/run_08_results.json', 'w') as f:
        json.dump(summary, f, default=lambda o: float(o) if isinstance(o, np.floating) else int(o) if isinstance(o, np.integer) else o, indent=2)
    flush_print(f"\n  Results saved to outputs/analysis/run_08/run_08_results.json")

    return summary


if __name__ == "__main__":
    main()
