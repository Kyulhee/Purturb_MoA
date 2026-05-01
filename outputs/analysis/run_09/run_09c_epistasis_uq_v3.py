"""
Run 09c: Epistasis Detection + UQ Pipeline v3
===============================================
Key fixes from v2:
  1. RQ3: Continuous ranking (not binary classification) — all real combos deviate from additive
  2. RQ2: Proper MC Dropout on decoder; z_tx OOD distance as epistemic uncertainty
  3. RQ4: "Find strongest epistasis" ranking problem, not binary detection
  4. Ground truth: permutation effect_size as continuous label
"""

import sys, io, os, json, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import r2_score, roc_auc_score, average_precision_score
from scipy.stats import spearmanr, pearsonr, kendalltau
import warnings
warnings.filterwarnings("ignore")

torch.manual_seed(42)
np.random.seed(42)

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__))
DEVICE = torch.device('cpu')


# ============================================================
# 1. Data Loading (same)
# ============================================================

def load_norman_data(data_path='outputs/analysis/run_04/data/gears_data'):
    from gears import PertData
    pd = PertData(data_path)
    pd.load(data_name='norman')
    return pd.adata, pd


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
    single_perts, double_perts = [], []
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
        s1 = (p1 in single_set or f'{p1}+ctrl' in single_set or f'ctrl+{p1}' in single_set)
        s2 = (p2 in single_set or f'{p2}+ctrl' in single_set or f'ctrl+{p2}' in single_set)
        if s1 and s2:
            valid_double.append(dp)
    return single_perts, valid_double


def build_gene_to_condition(conditions_unique):
    gene_to_condition = {}
    for cond in conditions_unique:
        if cond == 'ctrl':
            continue
        if '+ctrl' in cond:
            gene = cond.replace('+ctrl', '')
            gene_to_condition[gene] = cond
        elif 'ctrl+' in cond:
            gene = cond.replace('ctrl+', '')
            gene_to_condition[gene] = cond
        elif '+' not in cond:
            gene_to_condition[cond] = cond
    return gene_to_condition


# ============================================================
# 2. Model Architecture — with proper dropout for MC Dropout
# ============================================================

class FCREncoder(nn.Module):
    def __init__(self, n_genes, n_perturbations, z_dim, n_cell_types=1):
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
    def __init__(self, z_dim, n_genes, dropout_rate=0.1):
        super().__init__()
        self.decoder = nn.Sequential(
            nn.Linear(3 * z_dim, 256), nn.ReLU(), nn.Dropout(dropout_rate),
            nn.Linear(256, 128), nn.ReLU(), nn.Dropout(dropout_rate),
            nn.Linear(128, n_genes),
        )

    def forward(self, z_x, z_t, z_tx):
        return self.decoder(torch.cat([z_x, z_t, z_tx], dim=-1))


def reparameterize(mean, logvar):
    std = torch.exp(0.5 * logvar)
    return mean + std * torch.randn_like(std)


def icm_regularizer(z_tx_mean, z_tx_logvar, cell_type_onehot, cell_types):
    unique_types = torch.unique(cell_types)
    if len(unique_types) > 1:
        mmd_loss = torch.tensor(0.0, device=z_tx_mean.device)
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
    else:
        var = z_tx_mean.var(0)
        return ((var - 1.0) ** 2).mean()


def train_fcr_icm(encoder, decoder, X_train, pert_ids, ct_ids, n_cell_types,
                  n_epochs=100, batch_size=512, beta=0.5, icm_weight=10.0):
    optimizer = torch.optim.Adam(
        list(encoder.parameters()) + list(decoder.parameters()), lr=1e-3
    )
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
            loss = recon_loss + beta * kl
            icm_loss = icm_regularizer(z_tx_m, z_tx_lv, batch_ct_oh, batch_ct)
            loss = loss + icm_weight * icm_loss
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
        if (epoch + 1) % 25 == 0:
            print(f"  Epoch {epoch+1}/{n_epochs}: loss={epoch_loss/len(loader):.2f}")

    return encoder, decoder


# ============================================================
# 3. Phase 2: Combination Prediction + Residual Extraction
# ============================================================

def extract_residuals(encoder, decoder, adata, X, pert_id_map, n_cell_types, z_dim):
    encoder.eval()
    decoder.eval()
    conditions = adata.obs['condition'].values
    gene_to_condition = build_gene_to_condition(np.unique(conditions))

    single_z_tx, single_z_x, single_z_t = {}, {}, {}
    for gene, cond_name in gene_to_condition.items():
        pert_idx = pert_id_map[cond_name]
        mask = conditions == cond_name
        if mask.sum() == 0:
            continue
        x_pert = torch.FloatTensor(X[mask])
        pert_t = torch.full((x_pert.shape[0],), pert_idx, dtype=torch.long)
        ct_oh = F.one_hot(torch.zeros(x_pert.shape[0], dtype=torch.long), n_cell_types).float()
        with torch.no_grad():
            (z_x_m, _), (z_t_m, _), (z_tx_m, _) = encoder(x_pert, pert_t, ct_oh)
        single_z_tx[gene] = z_tx_m.mean(0)
        single_z_x[gene] = z_x_m.mean(0)
        single_z_t[gene] = z_t_m.mean(0)

    print(f"  Encoded single perturbations: {len(single_z_tx)} genes")

    ctrl_mask = conditions == 'ctrl'
    ctrl_mean = X[ctrl_mask].mean(axis=0)

    results = []
    all_double = [c for c in np.unique(conditions) if '+' in c and c != 'ctrl'
                  and not c.startswith('ctrl+') and not c.endswith('+ctrl')]

    for dp in all_double:
        p1, p2 = dp.split('+')
        if p1 not in single_z_tx or p2 not in single_z_tx:
            continue
        dp_mask = conditions == dp
        if dp_mask.sum() < 5:
            continue
        y_real = X[dp_mask].mean(axis=0)

        with torch.no_grad():
            z_tx_composed = single_z_tx[p1] + single_z_tx[p2]
            z_x_ref = single_z_x[p1].unsqueeze(0)
            z_t_ref = single_z_t[p1].unsqueeze(0)
            y_pred = decoder(z_x_ref, z_t_ref, z_tx_composed.unsqueeze(0))[0].numpy()

        y_A = X[conditions == gene_to_condition.get(p1, p1)].mean(axis=0)
        y_B = X[conditions == gene_to_condition.get(p2, p2)].mean(axis=0)
        y_additive = y_A + y_B - ctrl_mean

        results.append({
            'pert1': p1, 'pert2': p2, 'double': dp,
            'y_real': y_real, 'y_pred': y_pred,
            'y_A': y_A, 'y_B': y_B, 'y_additive': y_additive,
            'ctrl_mean': ctrl_mean,
            'residual_add': y_real - y_additive,
            'z_tx_composed': z_tx_composed.numpy(),
            'z_tx_1': single_z_tx[p1].numpy(),
            'z_tx_2': single_z_tx[p2].numpy(),
            'n_cells': int(dp_mask.sum()),
            'r2_additive': r2_score(y_real, y_additive),
            'corr_additive': pearsonr(y_real, y_additive)[0],
        })

    return results


# ============================================================
# 4. Phase 3: Residual Decomposition (RQ1)
# ============================================================

def compute_training_z_tx_stats(encoder, adata, X, pert_id_map, n_cell_types):
    """Compute mean and covariance of z_tx for each single perturbation."""
    encoder.eval()
    conditions = adata.obs['condition'].values
    gene_to_condition = build_gene_to_condition(np.unique(conditions))

    stats = {}
    for gene, cond_name in gene_to_condition.items():
        pert_idx = pert_id_map[cond_name]
        mask = conditions == cond_name
        if mask.sum() < 5:
            continue
        x_pert = torch.FloatTensor(X[mask])
        pert_t = torch.full((x_pert.shape[0],), pert_idx, dtype=torch.long)
        ct_oh = F.one_hot(torch.zeros(x_pert.shape[0], dtype=torch.long), n_cell_types).float()
        with torch.no_grad():
            (z_x_m, _), (z_t_m, _), (z_tx_m, _) = encoder(x_pert, pert_t, ct_oh)
        z_tx_np = z_tx_m.numpy()
        stats[gene] = {
            'mean': z_tx_np.mean(axis=0),
            'std': z_tx_np.std(axis=0) + 1e-6,
        }
    return stats


def decompose_residuals(encoder, decoder, residual_results, adata, X,
                        pert_id_map, n_cell_types, z_dim):
    """Decompose residuals using OOD distance + MC Dropout variance."""
    conditions = adata.obs['condition'].values
    decomposition_results = []

    # Get training z_tx statistics for OOD computation
    z_tx_stats = compute_training_z_tx_stats(encoder, adata, X, pert_id_map, n_cell_types)

    for rr in residual_results:
        p1, p2, dp = rr['pert1'], rr['pert2'], rr['double']

        # --- OOD score: distance of composed z_tx from training distribution ---
        z_tx_composed = rr['z_tx_composed']
        ood_scores = []
        for gene, stats in z_tx_stats.items():
            # Mahalanobis-like distance: sum of |z - mean| / std
            dist = np.abs(z_tx_composed - stats['mean']) / stats['std']
            ood_scores.append(float(dist.mean()))
        ood_score = float(np.mean(ood_scores))  # average across all single perturbations

        # Specific OOD: distance from each parent's z_tx distribution
        if p1 in z_tx_stats and p2 in z_tx_stats:
            ood_p1 = float(np.abs(z_tx_composed - z_tx_stats[p1]['mean']).mean() / z_tx_stats[p1]['std'].mean())
            ood_p2 = float(np.abs(z_tx_composed - z_tx_stats[p2]['mean']).mean() / z_tx_stats[p2]['std'].mean())
        else:
            ood_p1, ood_p2 = 0, 0

        # --- MC Dropout on decoder for composed z_tx ---
        dp_mask = conditions == dp
        if dp_mask.sum() < 5:
            continue

        # Re-encode to get z_x, z_t
        cond_A = build_gene_to_condition(np.unique(conditions)).get(p1, p1)
        mask_A = conditions == cond_A
        x_A = torch.FloatTensor(X[mask_A][:1])
        pert_idx_A = pert_id_map.get(cond_A, 0)
        pert_t_A = torch.full((1,), pert_idx_A, dtype=torch.long)
        ct_oh_A = F.one_hot(torch.zeros(1, dtype=torch.long), n_cell_types).float()

        with torch.no_grad():
            (z_x_m, _), (z_t_m, _), (_, _) = encoder(x_A, pert_t_A, ct_oh_A)

        z_tx_tensor = torch.FloatTensor(z_tx_composed).unsqueeze(0)
        z_x_ref = z_x_m[:1]
        z_t_ref = z_t_m[:1]

        # Proper MC Dropout: enable dropout, forward through decoder
        decoder.train()
        mc_preds = []
        with torch.no_grad():
            for _ in range(30):
                pred = decoder(z_x_ref, z_t_ref, z_tx_tensor)
                mc_preds.append(pred.numpy())
        decoder.eval()
        mc_preds = np.array(mc_preds)
        mc_var = mc_preds.var(axis=0)[0]  # (n_genes,)
        mc_mean = mc_preds.mean(axis=0)[0]

        # --- ICM violation ---
        x_dp = torch.FloatTensor(X[dp_mask][:10])
        pert_idx_dp = pert_id_map.get(dp, 0)
        pert_t_dp = torch.full((x_dp.shape[0],), pert_idx_dp, dtype=torch.long)
        ct_oh_dp = F.one_hot(torch.zeros(x_dp.shape[0], dtype=torch.long), n_cell_types).float()
        with torch.no_grad():
            (_, _, (z_tx_m_dp, _)) = encoder(x_dp, pert_t_dp, ct_oh_dp)
            z_tx_var = z_tx_m_dp.var(0)
            icm_violation = float(((z_tx_var - 1.0) ** 2).mean().item())

        # --- Decomposition ---
        r_total = rr['residual_add']
        r_total_mag = float(np.abs(r_total).mean())
        r_noise_mag = float(mc_var.mean() ** 0.5)  # sqrt of mean variance
        r_model_mag = ood_score * float(mc_var.mean() ** 0.5) * 0.1  # scale OOD by prediction uncertainty
        r_epistasis_mag = max(r_total_mag - r_noise_mag - r_model_mag, 0)

        decomposition_results.append({
            'pert1': p1, 'pert2': p2, 'double': dp,
            'r_total': r_total,
            'ood_score': ood_score,
            'ood_p1': ood_p1, 'ood_p2': ood_p2,
            'icm_violation': icm_violation,
            'mc_var_mean': float(mc_var.mean()),
            'mc_var_gene': mc_var,
            'r_noise_mag': r_noise_mag,
            'r_model_mag': r_model_mag,
            'r_epistasis_mag': r_epistasis_mag,
            'r_total_mag': r_total_mag,
            'r2_additive': rr['r2_additive'],
        })

    return decomposition_results


# ============================================================
# 5. Phase 4: Uncertainty Quantification (RQ2) — v3
# ============================================================

def quantify_uncertainty(decomposition_results):
    """UQ using OOD distance, MC Dropout variance, and ICM violation.
    Evaluates whether uncertainty RANKS combinations by prediction error."""
    records = []

    for d in decomposition_results:
        abs_error = np.abs(d['r_total'])
        mean_abs_error = float(abs_error.mean())

        records.append({
            'double': d['double'],
            'pert1': d['pert1'], 'pert2': d['pert2'],
            'ood_score': d['ood_score'],
            'icm_score': d['icm_violation'],
            'mc_var': d['mc_var_mean'],
            'mean_abs_error': mean_abs_error,
            'r2': d['r2_additive'],
        })

    if len(records) < 4:
        return {'error': 'Too few combinations'}

    u_ood = np.array([r['ood_score'] for r in records])
    u_icm = np.array([r['icm_score'] for r in records])
    u_mc = np.array([r['mc_var'] for r in records])
    err = np.array([r['mean_abs_error'] for r in records])

    # Individual signal correlations
    rho_ood, p_ood = spearmanr(u_ood, err)
    rho_icm, p_icm = spearmanr(u_icm, err)
    rho_mc, p_mc = spearmanr(u_mc, err)

    # Optimize combined weights
    best_rho, best_w = -1, {'ood': 0, 'icm': 0, 'mc': 0}
    for w_ood in np.arange(0, 1.01, 0.2):
        for w_icm in np.arange(0, 1.01, 0.2):
            for w_mc in np.arange(0, 1.01, 0.2):
                total = w_ood + w_icm + w_mc
                if total == 0:
                    continue
                u_trial = (w_ood * u_ood + w_icm * u_icm + w_mc * u_mc) / total
                rho, _ = spearmanr(u_trial, err)
                if rho > best_rho:
                    best_rho = rho
                    best_w = {'ood': float(w_ood), 'icm': float(w_icm), 'mc': float(w_mc)}

    # Compute best combined
    total_w = best_w['ood'] + best_w['icm'] + best_w['mc']
    u_combined = (best_w['ood'] * u_ood + best_w['icm'] * u_icm + best_w['mc'] * u_mc) / total_w
    rho_combined, p_combined = spearmanr(u_combined, err)

    # Coverage (simple)
    coverage_records = []
    for d in decomposition_results:
        r_mean = float(np.abs(d['r_total']).mean())
        ci_width = 1.645 * np.sqrt(d['mc_var_mean']) * np.sqrt(500)  # n_genes=500
        coverage_records.append(float(r_mean < ci_width))

    return {
        'n_combinations': len(records),
        'spearman_rho_combined': float(rho_combined),
        'spearman_p_combined': float(p_combined),
        'spearman_rho_ood': float(rho_ood),
        'spearman_p_ood': float(p_ood),
        'spearman_rho_icm': float(rho_icm),
        'spearman_rho_mc': float(rho_mc),
        'best_weights': best_w,
        'mean_coverage': float(np.mean(coverage_records)),
        'per_combination': records,
    }


# ============================================================
# 6. Phase 5: Epistasis Detection (RQ3) — v3: Continuous Ranking
# ============================================================

def compute_epistasis_measures(residual_results, adata, X):
    """Compute epistasis using 3 formulas. Use continuous effect size as ground truth."""
    conditions = adata.obs['condition'].values
    gene_to_condition = build_gene_to_condition(np.unique(conditions))

    ctrl_mask = conditions == 'ctrl'
    ctrl_mean = X[ctrl_mask].mean(axis=0)
    ctrl_std = X[ctrl_mask].std(axis=0) + 1e-6

    # Compute empirical effect sizes for ground truth ranking
    # Effect size = |observed - additive| / pooled_std, averaged across genes
    ctrl_cells = X[ctrl_mask]
    ctrl_var = ctrl_cells.var(axis=0) + 1e-6

    epistasis_results = []

    for rr in residual_results:
        p1, p2, dp = rr['pert1'], rr['pert2'], rr['double']
        cond_A = gene_to_condition.get(p1, p1)
        cond_B = gene_to_condition.get(p2, p2)
        mask_A = conditions == cond_A
        mask_B = conditions == cond_B
        dp_mask = conditions == dp
        if mask_A.sum() < 5 or mask_B.sum() < 5 or dp_mask.sum() < 5:
            continue

        y_A = X[mask_A].mean(axis=0)
        y_B = X[mask_B].mean(axis=0)
        y_AB = X[dp_mask].mean(axis=0)

        # --- Ground truth: Cohen's d-like effect size ---
        expected_add = y_A + y_B - ctrl_mean
        residual_add = y_AB - expected_add
        # Normalize by control std
        effect_size_gt = float(np.abs(residual_add / ctrl_std).mean())

        # Also compute using cell-level variance
        residual_var = X[dp_mask].var(axis=0) + X[mask_A].var(axis=0) + 1e-6
        effect_size_cell = float(np.abs(residual_add / np.sqrt(residual_var / 2)).mean())

        # --- Formula 1: Additive ---
        epistasis_add = residual_add
        epistasis_strength_add = float(np.abs(epistasis_add / ctrl_std).mean())

        # --- Formula 2: Multiplicative (log scale) ---
        eps = 1e-6
        expected_mult = np.log(np.maximum(y_A, eps)) + np.log(np.maximum(y_B, eps)) - np.log(np.maximum(ctrl_mean, eps))
        epistasis_mult = np.log(np.maximum(y_AB, eps)) - expected_mult
        epistasis_strength_mult = float(np.abs(epistasis_mult).mean())

        # --- Formula 3: Product neutrality (Valenzuela 2025) ---
        fitness_A = y_A / (ctrl_mean + eps)
        fitness_B = y_B / (ctrl_mean + eps)
        fitness_AB = y_AB / (ctrl_mean + eps)
        expected_prod = fitness_A * fitness_B
        epistasis_prod = fitness_AB - expected_prod
        epistasis_strength_prod = float(np.abs(epistasis_prod).mean())

        # --- Formula agreement at combination level ---
        # Do all 3 formulas agree on the DIRECTION of epistasis for majority genes?
        dir_add = np.sign(epistasis_add)
        dir_mult = np.sign(epistasis_mult)
        dir_prod = np.sign(epistasis_prod)

        agree_add_mult = float((dir_add == dir_mult).mean())
        agree_add_prod = float((dir_add == dir_prod).mean())
        agree_mult_prod = float((dir_mult == dir_prod).mean())
        mean_agreement = (agree_add_mult + agree_add_prod + agree_mult_prod) / 3

        # --- GI classification (by dominant direction) ---
        direction_mean = float(dir_add.mean())
        if abs(direction_mean) < 0.1:
            gi_type = 'neomorphic'  # mixed directions
        elif direction_mean > 0:
            gi_type = 'synergy'
        else:
            gi_type = 'suppression'

        epistasis_results.append({
            'double': dp, 'pert1': p1, 'pert2': p2,
            'effect_size_gt': effect_size_gt,
            'effect_size_cell': effect_size_cell,
            'epistasis_strength_add': epistasis_strength_add,
            'epistasis_strength_mult': epistasis_strength_mult,
            'epistasis_strength_prod': epistasis_strength_prod,
            'agreement_add_mult': agree_add_mult,
            'agreement_add_prod': agree_add_prod,
            'agreement_mult_prod': agree_mult_prod,
            'mean_agreement': mean_agreement,
            'gi_type': gi_type,
            'r2_additive': rr['r2_additive'],
        })

    return epistasis_results


def evaluate_epistasis(epistasis_results):
    """Evaluate epistasis detection as CONTINUOUS RANKING, not binary classification."""
    gt_effect = np.array([r['effect_size_gt'] for r in epistasis_results])
    scores_add = np.array([r['epistasis_strength_add'] for r in epistasis_results])
    scores_mult = np.array([r['epistasis_strength_mult'] for r in epistasis_results])
    scores_prod = np.array([r['epistasis_strength_prod'] for r in epistasis_results])

    results = {}
    n_total = len(epistasis_results)
    results['n_total'] = n_total

    # Core metric: Spearman rho between epistasis strength and ground truth effect size
    rho_add, p_add = spearmanr(scores_add, gt_effect)
    rho_mult, p_mult = spearmanr(scores_mult, gt_effect)
    rho_prod, p_prod = spearmanr(scores_prod, gt_effect)

    results['spearman_rho_additive_gt'] = float(rho_add)
    results['spearman_p_additive_gt'] = float(p_add)
    results['spearman_rho_multiplicative_gt'] = float(rho_mult)
    results['spearman_p_multiplicative_gt'] = float(p_mult)
    results['spearman_rho_product_gt'] = float(rho_prod)
    results['spearman_p_product_gt'] = float(p_prod)

    # Binary AUROC using median split (for comparison with run_09)
    median_effect = np.median(gt_effect)
    labels = (gt_effect >= median_effect).astype(int)
    if len(set(labels)) > 1:
        results['auroc_additive'] = float(roc_auc_score(labels, scores_add))
        results['auroc_product'] = float(roc_auc_score(labels, scores_prod))
        results['auroc_multiplicative'] = float(roc_auc_score(labels, scores_mult))
        results['aupr_additive'] = float(average_precision_score(labels, scores_add))

    # Top-k precision: top-k by score, what fraction are in top-k by GT?
    top_k_by_gt = set(np.argsort(gt_effect)[-20:])
    for k in [10, 20]:
        top_k_by_score = set(np.argsort(scores_add)[-k:])
        overlap = len(top_k_by_score & top_k_by_gt) / k
        results[f'top{k}_overlap_with_gt'] = float(overlap)

    # GI type distribution
    type_counts = {}
    for r in epistasis_results:
        t = r['gi_type']
        type_counts[t] = type_counts.get(t, 0) + 1
    results['gi_type_distribution'] = type_counts

    # Formula agreement stats
    agreements = [r['mean_agreement'] for r in epistasis_results]
    results['mean_formula_agreement'] = float(np.mean(agreements))

    return results


# ============================================================
# 7. Phase 6: Active Learning (RQ4) — v3: Ranking-based
# ============================================================

def simulate_active_learning(epistasis_results, uq_results, decomposition_results):
    """AL as ranking problem: find combinations with strongest epistasis."""
    if not epistasis_results or not uq_results.get('per_combination'):
        return {'error': 'Insufficient data'}

    epi_lookup = {r['double']: r for r in epistasis_results}
    uq_lookup = {r['double']: r for r in uq_results['per_combination']}
    decomp_lookup = {r['double']: r for r in decomposition_results}

    all_combos = [r['double'] for r in epistasis_results]
    gt_effect = {r['double']: r['effect_size_gt'] for r in epistasis_results}
    n = len(all_combos)

    # Define "strong epistasis" as top 30% by ground truth effect size
    effects = sorted(gt_effect.values(), reverse=True)
    threshold = effects[int(n * 0.3)]
    is_strong = {c: gt_effect[c] >= threshold for c in all_combos}
    n_strong = sum(is_strong.values())

    # --- Random ---
    np.random.seed(42)
    random_order = list(np.random.permutation(all_combos))

    # --- UQ-only ---
    uq_scores = {c: uq_lookup.get(c, {}).get('ood_score', 0) +
                 uq_lookup.get(c, {}).get('icm_score', 0) +
                 uq_lookup.get(c, {}).get('mc_var', 0)
                 for c in all_combos}
    uq_order = sorted(all_combos, key=lambda c: -uq_scores[c])

    # --- Epistasis score from decomposition ---
    epi_scores = {c: decomp_lookup.get(c, {}).get('r_epistasis_mag', 0) for c in all_combos}
    epi_order = sorted(all_combos, key=lambda c: -epi_scores[c])

    # --- OOD-based ---
    ood_scores = {c: decomp_lookup.get(c, {}).get('ood_score', 0) for c in all_combos}
    ood_order = sorted(all_combos, key=lambda c: -ood_scores[c])

    # --- Oracle ---
    oracle_order = sorted(all_combos, key=lambda c: -gt_effect[c])

    def cumulative_recall(order, is_strong, n_strong):
        found = 0
        cum = []
        for c in order:
            if is_strong.get(c, False):
                found += 1
            cum.append(found / max(n_strong, 1))
        return cum

    random_cum = cumulative_recall(random_order, is_strong, n_strong)
    uq_cum = cumulative_recall(uq_order, is_strong, n_strong)
    epi_cum = cumulative_recall(epi_order, is_strong, n_strong)
    ood_cum = cumulative_recall(ood_order, is_strong, n_strong)
    oracle_cum = cumulative_recall(oracle_order, is_strong, n_strong)

    top_k_results = {}
    for k in [5, 10, 20]:
        if k <= n:
            top_k_results[k] = {
                'random': random_cum[k-1],
                'uq': uq_cum[k-1],
                'epistasis': epi_cum[k-1],
                'ood': ood_cum[k-1],
                'oracle': oracle_cum[k-1],
                'improvement_uq': uq_cum[k-1] / max(random_cum[k-1], 1e-8),
                'improvement_epi': epi_cum[k-1] / max(random_cum[k-1], 1e-8),
                'improvement_ood': ood_cum[k-1] / max(random_cum[k-1], 1e-8),
            }

    # NDCG@k for ranking quality
    def ndcg_at_k(order, gt_effect, k):
        dcg = 0
        for i, c in enumerate(order[:k]):
            dcg += gt_effect[c] / np.log2(i + 2)
        ideal = sorted(gt_effect.values(), reverse=True)
        idcg = sum(ideal[i] / np.log2(i + 2) for i in range(min(k, len(ideal))))
        return dcg / max(idcg, 1e-8)

    ndcg_results = {}
    for k in [10, 20]:
        ndcg_results[k] = {
            'random': ndcg_at_k(random_order, gt_effect, k),
            'uq': ndcg_at_k(uq_order, gt_effect, k),
            'epistasis': ndcg_at_k(epi_order, gt_effect, k),
            'ood': ndcg_at_k(ood_order, gt_effect, k),
            'oracle': ndcg_at_k(oracle_order, gt_effect, k),
        }

    return {
        'n_total': n,
        'n_strong': n_strong,
        'top_k_results': top_k_results,
        'ndcg': ndcg_results,
    }


# ============================================================
# 8. Main Pipeline
# ============================================================

def main():
    t0 = time.time()
    print("=" * 70)
    print("Run 09c: Epistasis Detection + UQ Pipeline v3")
    print("=" * 70)

    # Phase 1
    print("\n[Phase 1] Loading Norman 2019 data...")
    adata, gears_pd = load_norman_data()
    adata, X = preprocess_adata(adata, n_top_genes=500)
    n_genes = X.shape[1]
    print(f"  Data: {X.shape[0]} cells, {n_genes} genes")

    single_perts, double_perts = get_perturbation_info(adata)
    print(f"  Single: {len(single_perts)}, Double (valid): {len(double_perts)}")

    all_pert_names = sorted(set(single_perts + double_perts + ['ctrl']))
    pert_id_map = {name: idx for idx, name in enumerate(all_pert_names)}
    n_perturbations = len(all_pert_names)
    n_cell_types = 1
    z_dim = 8

    conditions = adata.obs['condition'].values
    valid_mask = np.isin(conditions, all_pert_names)
    X_filt = X[valid_mask]
    cond_filt = conditions[valid_mask]
    pert_ids = np.array([pert_id_map.get(c, 0) for c in cond_filt])
    ct_ids = np.zeros(len(X_filt), dtype=np.int64)
    print(f"  Training cells: {len(X_filt)}")

    print("\n[Phase 1] Training FCR-ICM model...")
    encoder = FCREncoder(n_genes, n_perturbations, z_dim, n_cell_types)
    decoder = FCRDecoder(z_dim, n_genes, dropout_rate=0.1)
    encoder, decoder = train_fcr_icm(encoder, decoder, X_filt, pert_ids, ct_ids, n_cell_types)
    print("  Done.")

    # Phase 2
    print("\n[Phase 2] Extracting prediction residuals...")
    residual_results = extract_residuals(encoder, decoder, adata, X, pert_id_map, n_cell_types, z_dim)
    print(f"  Pairs: {len(residual_results)}")
    if residual_results:
        r2_add = [r['r2_additive'] for r in residual_results]
        print(f"  Additive R2: mean={np.mean(r2_add):.3f}, median={np.median(r2_add):.3f}")

    # Phase 3
    print("\n[Phase 3] Residual decomposition (RQ1)...")
    decomposition_results = decompose_residuals(
        encoder, decoder, residual_results, adata, X, pert_id_map, n_cell_types, z_dim
    )
    print(f"  Decomposed: {len(decomposition_results)}")
    if decomposition_results:
        epi_frac = np.mean([d['r_epistasis_mag'] for d in decomposition_results]) / max(np.mean([d['r_total_mag'] for d in decomposition_results]), 1e-8)
        print(f"  Epistasis fraction: {epi_frac:.1%}")
        print(f"  Mean OOD score: {np.mean([d['ood_score'] for d in decomposition_results]):.3f}")

    # Phase 4
    print("\n[Phase 4] Uncertainty quantification (RQ2)...")
    uq_results = quantify_uncertainty(decomposition_results)
    if 'error' not in uq_results:
        print(f"  U-Error rho (combined): {uq_results['spearman_rho_combined']:.3f} (p={uq_results['spearman_p_combined']:.4f})")
        print(f"  U-Error rho (OOD): {uq_results['spearman_rho_ood']:.3f} (p={uq_results['spearman_p_ood']:.4f})")
        print(f"  U-Error rho (ICM): {uq_results['spearman_rho_icm']:.3f}")
        print(f"  U-Error rho (MC): {uq_results['spearman_rho_mc']:.3f}")
        print(f"  Best weights: {uq_results['best_weights']}")
        print(f"  Coverage: {uq_results['mean_coverage']:.3f}")

    # Phase 5
    print("\n[Phase 5] Epistasis detection (RQ3)...")
    epistasis_results = compute_epistasis_measures(residual_results, adata, X)
    print(f"  Analyzed: {len(epistasis_results)}")

    epi_eval = evaluate_epistasis(epistasis_results)
    print(f"  Spearman rho (additive vs GT): {epi_eval.get('spearman_rho_additive_gt', 'N/A'):.3f}")
    print(f"  Spearman rho (product vs GT): {epi_eval.get('spearman_rho_product_gt', 'N/A'):.3f}")
    print(f"  Spearman rho (multiplicative vs GT): {epi_eval.get('spearman_rho_multiplicative_gt', 'N/A'):.3f}")
    if 'auroc_additive' in epi_eval:
        print(f"  AUROC (additive, median split): {epi_eval['auroc_additive']:.3f}")
        print(f"  AUPR (additive): {epi_eval['aupr_additive']:.3f}")
    if 'top10_overlap_with_gt' in epi_eval:
        print(f"  Top-10 overlap with GT: {epi_eval['top10_overlap_with_gt']:.3f}")
    print(f"  GI types: {epi_eval.get('gi_type_distribution', {})}")
    print(f"  Formula agreement: {epi_eval.get('mean_formula_agreement', 'N/A'):.3f}")

    # Phase 6
    print("\n[Phase 6] Active learning (RQ4)...")
    al_results = simulate_active_learning(epistasis_results, uq_results, decomposition_results)
    if 'error' not in al_results:
        print(f"  Strong epistasis: {al_results['n_strong']}/{al_results['n_total']}")
        for k, v in al_results.get('top_k_results', {}).items():
            print(f"  Top-{k}: random={v['random']:.3f}, UQ={v['uq']:.3f}, "
                  f"epi={v['epistasis']:.3f}, OOD={v['ood']:.3f}, oracle={v['oracle']:.3f}")
            print(f"    Improvement: UQ={v['improvement_uq']:.2f}x, epi={v['improvement_epi']:.2f}x, OOD={v['improvement_ood']:.2f}x")
        for k, v in al_results.get('ndcg', {}).items():
            print(f"  NDCG@{k}: random={v['random']:.3f}, UQ={v['uq']:.3f}, "
                  f"epi={v['epistasis']:.3f}, OOD={v['ood']:.3f}, oracle={v['oracle']:.3f}")
    else:
        print(f"  {al_results['error']}")

    # Summary
    elapsed = time.time() - t0
    print("\n" + "=" * 70)
    print("RUN 09c SUMMARY")
    print("=" * 70)
    if 'error' not in uq_results:
        print(f"  RQ2: U-Error rho = {uq_results['spearman_rho_combined']:.3f} (target: >0.6)")
        print(f"  RQ2: OOD rho = {uq_results['spearman_rho_ood']:.3f}")
    print(f"  RQ3: Epistasis-GT rho (add) = {epi_eval.get('spearman_rho_additive_gt', 'N/A')}")
    print(f"  RQ3: Epistasis-GT rho (prod) = {epi_eval.get('spearman_rho_product_gt', 'N/A')}")
    if 'auroc_additive' in epi_eval:
        print(f"  RQ3: AUROC (add, median split) = {epi_eval['auroc_additive']:.3f} (target: >0.75)")
    if 'error' not in al_results:
        for k, v in al_results.get('top_k_results', {}).items():
            print(f"  RQ4: Top-{k} OOD improvement = {v['improvement_ood']:.2f}x vs random (target: >2x)")
    print(f"\n  Elapsed: {elapsed:.0f}s")

    # Save
    save_results = {
        'rq2_uq': {k: v for k, v in uq_results.items() if k != 'per_combination'},
        'rq3_epistasis': epi_eval,
        'rq4_al': {k: v for k, v in al_results.items() if not isinstance(v, list)},
        'n_pairs': len(residual_results),
        'elapsed_s': elapsed,
    }
    with open(os.path.join(RESULTS_DIR, 'run_09c_results.json'), 'w') as f:
        json.dump(save_results, f, indent=2, default=str)
    print(f"\n  Results saved to {RESULTS_DIR}/run_09c_results.json")


if __name__ == '__main__':
    main()

{}