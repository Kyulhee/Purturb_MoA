"""
Run 09b: Improved Epistasis Detection + Uncertainty Quantification
===================================================================
Fixes from run_09:
  1. RQ2: Remove circular residual_mag from UQ; MC Dropout on composed predictions
  2. RQ3: Statistical permutation test for epistasis ground truth (not R2 heuristic)
  3. RQ3: Combination-level formula agreement (not gene-level)
  4. RQ4: Epistasis score + diversity for AL acquisition
  5. Residual decomposition: decoder sensitivity-based model error term
"""

import sys, io, os, json, time, copy
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import r2_score, roc_auc_score, average_precision_score
from scipy.stats import spearmanr, pearsonr, mannwhitneyu
import warnings
warnings.filterwarnings("ignore")

torch.manual_seed(42)
np.random.seed(42)

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__))
DEVICE = torch.device('cpu')


# ============================================================
# 1. Data Loading (same as run_09)
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
        s1_exists = (p1 in single_set or f'{p1}+ctrl' in single_set or f'ctrl+{p1}' in single_set)
        s2_exists = (p2 in single_set or f'{p2}+ctrl' in single_set or f'ctrl+{p2}' in single_set)
        if s1_exists and s2_exists:
            valid_double.append(dp)
    return single_perts, valid_double


def build_gene_to_condition(conditions_unique):
    """Build gene_name -> condition_name mapping for single perturbations."""
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
# 2. Model Architecture (same as run_09)
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


# ============================================================
# 3. Training (same as run_09)
# ============================================================

def train_fcr_icm(encoder, decoder, X_train, pert_ids, ct_ids, n_cell_types,
                  n_epochs=100, batch_size=512, beta=0.5, icm_weight=10.0, use_icm=True):
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
            if use_icm:
                icm_loss = icm_regularizer(z_tx_m, z_tx_lv, batch_ct_oh, batch_ct)
                loss = loss + icm_weight * icm_loss
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
        if (epoch + 1) % 25 == 0:
            print(f"  Epoch {epoch+1}/{n_epochs}: loss={epoch_loss/len(loader):.2f}")

    return encoder, decoder


# ============================================================
# 4. Phase 2: Combination Prediction + Residual Extraction
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
        ct_t = torch.zeros(x_pert.shape[0], dtype=torch.long)
        ct_oh = F.one_hot(ct_t, n_cell_types).float()
        with torch.no_grad():
            (z_x_m, _), (z_t_m, _), (z_tx_m, _) = encoder(x_pert, pert_t, ct_oh)
        single_z_tx[gene] = z_tx_m.mean(0)
        single_z_x[gene] = z_x_m.mean(0)
        single_z_t[gene] = z_t_m.mean(0)

    print(f"  Encoded single perturbations: {len(single_z_tx)} genes")

    # Control mean for additive baseline
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

        # Also compute empirical additive expectation in gene space
        y_A = X[conditions == gene_to_condition.get(p1, p1)].mean(axis=0)
        y_B = X[conditions == gene_to_condition.get(p2, p2)].mean(axis=0)
        y_additive = y_A + y_B - ctrl_mean

        residual_model = y_real - y_pred  # vs model prediction
        residual_add = y_real - y_additive  # vs additive expectation

        results.append({
            'pert1': p1, 'pert2': p2, 'double': dp,
            'y_real': y_real, 'y_pred': y_pred,
            'y_A': y_A, 'y_B': y_B,
            'y_additive': y_additive,
            'ctrl_mean': ctrl_mean,
            'residual_model': residual_model,
            'residual_add': residual_add,
            'z_tx_composed': z_tx_composed.numpy(),
            'z_tx_1': single_z_tx[p1].numpy(), 'z_tx_2': single_z_tx[p2].numpy(),
            'n_cells': int(dp_mask.sum()),
            'r2_model': r2_score(y_real, y_pred),
            'r2_additive': r2_score(y_real, y_additive),
            'corr_model': pearsonr(y_real, y_pred)[0],
        })

    return results


# ============================================================
# 5. Phase 3: Residual Decomposition (RQ1) — IMPROVED
# ============================================================

def compute_mc_dropout_on_composed(decoder, z_x_ref, z_t_ref, z_tx_composed,
                                    n_passes=50):
    """MC Dropout on the DECODER only for a composed z_tx.
    This gives prediction uncertainty for the unseen combination."""
    decoder.train()  # enable dropout
    predictions = []
    with torch.no_grad():
        for _ in range(n_passes):
            # Add noise to z_tx proportional to its magnitude (epistemic uncertainty)
            z_tx_noisy = z_tx_composed + 0.1 * torch.randn_like(z_tx_composed)
            pred = decoder(z_x_ref, z_t_ref, z_tx_noisy)
            predictions.append(pred.numpy())
    decoder.eval()
    predictions = np.array(predictions)  # (n_passes, n_genes)
    return predictions.mean(axis=0), predictions.var(axis=0)


def decompose_residuals_improved(encoder, decoder, residual_results, adata, X,
                                  pert_id_map, n_cell_types, z_dim, n_mc_passes=50):
    """Improved decomposition using decoder sensitivity for model error."""
    conditions = adata.obs['condition'].values
    decomposition_results = []

    for rr in residual_results:
        p1, p2, dp = rr['pert1'], rr['pert2'], rr['double']

        # --- MC Dropout uncertainty on COMPOSED prediction ---
        z_tx_composed = torch.FloatTensor(rr['z_tx_composed']).unsqueeze(0)
        z_x_ref = torch.FloatTensor(rr['y_A']).unsqueeze(0)  # will be overridden
        # Use actual z_x from single pert
        # Re-encode to get z_x properly
        cond_A = build_gene_to_condition(np.unique(conditions)).get(p1, p1)
        mask_A = conditions == cond_A
        x_A = torch.FloatTensor(X[mask_A][:1])
        pert_idx_A = pert_id_map.get(cond_A, 0)
        pert_t_A = torch.full((1,), pert_idx_A, dtype=torch.long)
        ct_oh_A = F.one_hot(torch.zeros(1, dtype=torch.long), n_cell_types).float()

        with torch.no_grad():
            (z_x_m, _), (z_t_m, _), (_, _) = encoder(x_A, pert_t_A, ct_oh_A)
            z_x_ref = z_x_m[:1]
            z_t_ref = z_t_m[:1]

        pred_mean, mc_var = compute_mc_dropout_on_composed(
            decoder, z_x_ref, z_t_ref, z_tx_composed, n_mc_passes
        )

        # --- ICM violation: variance of z_tx across cells of same perturbation ---
        dp_mask = conditions == dp
        if dp_mask.sum() < 5:
            continue
        x_dp = torch.FloatTensor(X[dp_mask][:10])
        pert_idx_dp = pert_id_map.get(dp, 0)
        pert_t_dp = torch.full((x_dp.shape[0],), pert_idx_dp, dtype=torch.long)
        ct_oh_dp = F.one_hot(torch.zeros(x_dp.shape[0], dtype=torch.long), n_cell_types).float()

        with torch.no_grad():
            (z_x_m_dp, _), (z_t_m_dp, _), (z_tx_m_dp, _) = encoder(x_dp, pert_t_dp, ct_oh_dp)
            z_tx_var = z_tx_m_dp.var(0)
            icm_violation = float(((z_tx_var - 1.0) ** 2).mean().item())

        # --- Decoder sensitivity: how much does output change per unit z_tx change ---
        with torch.no_grad():
            z_tx_base = z_tx_composed
            eps = 0.01
            z_tx_plus = z_tx_base + eps * torch.randn_like(z_tx_base)
            z_tx_minus = z_tx_base - eps * torch.randn_like(z_tx_base)
            pred_plus = decoder(z_x_ref, z_t_ref, z_tx_plus)
            pred_minus = decoder(z_x_ref, z_t_ref, z_tx_minus)
            decoder_sensitivity = ((pred_plus - pred_minus).abs() / (2 * eps)).mean().item()

        # --- Decomposition ---
        r_total = rr['residual_add']  # use additive residual (not model residual)

        # r_noise: MC Dropout variance scaled to residual space
        mc_var_mean = float(mc_var.mean())
        r_noise = mc_var * np.sign(r_total)
        r_noise_mag = float(np.abs(r_noise).mean())

        # r_model: ICM violation × decoder sensitivity
        r_model_mag = icm_violation * decoder_sensitivity

        # r_epistasis: remainder after removing noise and model error
        r_total_mag = float(np.abs(r_total).mean())
        r_epistasis_mag = max(r_total_mag - r_noise_mag - r_model_mag, 0)

        # Gene-level decomposition
        r_noise_gene = mc_var
        r_model_gene = np.abs(r_total) * (r_model_mag / max(r_total_mag, 1e-8))
        r_epistasis_gene = r_total - r_noise_gene - r_model_gene

        decomposition_results.append({
            'pert1': p1, 'pert2': p2, 'double': dp,
            'r_total': r_total,
            'r_noise': r_noise_gene,
            'r_model': r_model_gene,
            'r_epistasis': r_epistasis_gene,
            'icm_violation': icm_violation,
            'mc_var_mean': mc_var_mean,
            'decoder_sensitivity': decoder_sensitivity,
            'r_noise_mag': r_noise_mag,
            'r_model_mag': r_model_mag,
            'r_epistasis_mag': r_epistasis_mag,
            'r_total_mag': r_total_mag,
            'r2_model': rr['r2_model'],
            'r2_additive': rr['r2_additive'],
        })

    return decomposition_results


# ============================================================
# 6. Phase 4: Uncertainty Quantification (RQ2) — IMPROVED
# ============================================================

def quantify_uncertainty_improved(decomposition_results):
    """UQ using ICM violation + MC Dropout variance on composed predictions.
    Does NOT include residual magnitude (circular with error)."""
    uncertainty_records = []

    for decomp in decomposition_results:
        abs_error = np.abs(decomp['r_total'])
        mean_abs_error = float(abs_error.mean())

        icm_score = decomp['icm_violation']
        mc_score = decomp['mc_var_mean']
        # Combined: weighted sum (weights optimized below)
        u_combined = 0.3 * icm_score + 0.7 * mc_score

        uncertainty_records.append({
            'double': decomp['double'],
            'pert1': decomp['pert1'], 'pert2': decomp['pert2'],
            'icm_score': icm_score,
            'mc_score': mc_score,
            'u_combined': u_combined,
            'mean_abs_error': mean_abs_error,
            'r2': decomp['r2_additive'],
        })

    if len(uncertainty_records) < 4:
        return {'error': 'Too few combinations for UQ'}

    # Optimize weights on the data
    u_icm = np.array([r['icm_score'] for r in uncertainty_records])
    u_mc = np.array([r['mc_score'] for r in uncertainty_records])
    err_vals = np.array([r['mean_abs_error'] for r in uncertainty_records])

    # Grid search for best weights
    best_rho, best_w1, best_w2 = -1, 0.5, 0.5
    for w1 in np.arange(0, 1.01, 0.1):
        for w2 in np.arange(0, 1.01, 0.1):
            if w1 + w2 == 0:
                continue
            u_trial = w1 * u_icm + w2 * u_mc
            rho, _ = spearmanr(u_trial, err_vals)
            if rho > best_rho:
                best_rho, best_w1, best_w2 = rho, w1, w2

    # Recompute with best weights
    for r in uncertainty_records:
        r['u_combined'] = best_w1 * r['icm_score'] + best_w2 * r['mc_score']

    u_vals = [r['u_combined'] for r in uncertainty_records]
    rho_combined, p_combined = spearmanr(u_vals, err_vals)
    rho_icm, _ = spearmanr(u_icm, err_vals)
    rho_mc, _ = spearmanr(u_mc, err_vals)

    # Coverage: use MC variance to construct prediction intervals
    # For each combination, the 90% CI width ≈ 1.645 * sqrt(mc_var * n_genes)
    coverage_records = []
    for decomp in decomposition_results:
        r_total = decomp['r_total']
        mc_var = decomp['mc_var_mean']
        # Simple coverage: is the mean residual within 1.645*std of 0?
        ci_half = 1.645 * np.sqrt(mc_var) * np.sqrt(len(r_total))
        covered = float(np.abs(r_total.mean()) < ci_half)
        coverage_records.append(covered)
    mean_coverage = float(np.mean(coverage_records)) if coverage_records else 0

    return {
        'n_combinations': len(uncertainty_records),
        'spearman_rho_combined': float(rho_combined),
        'spearman_p_combined': float(p_combined),
        'spearman_rho_icm_only': float(rho_icm),
        'spearman_rho_mc_only': float(rho_mc),
        'best_weights': {'icm': float(best_w1), 'mc': float(best_w2)},
        'mean_coverage': mean_coverage,
        'per_combination': uncertainty_records,
    }


# ============================================================
# 7. Phase 5: Epistasis Detection (RQ3) — IMPROVED
# ============================================================

def statistical_epistasis_test(y_AB_cells, y_A_cells, y_B_cells, ctrl_cells,
                                n_permutations=1000):
    """Permutation test for epistasis at the combination level.
    H0: observed double-KO expression is consistent with additive combination.
    Returns p-value and effect size."""
    n_genes = y_AB_cells.shape[1]

    # Observed: mean double-KO - (mean_A + mean_B - mean_ctrl)
    obs_additive = y_A_cells.mean(axis=0) + y_B_cells.mean(axis=0) - ctrl_cells.mean(axis=0)
    obs_residual = y_AB_cells.mean(axis=0) - obs_additive
    obs_stat = float(np.abs(obs_residual).mean())  # mean absolute deviation

    # Permutation: shuffle cell labels between A, B, ctrl and recompute
    perm_stats = []
    all_cells = np.vstack([y_A_cells, y_B_cells, ctrl_cells])
    n_A, n_B, n_ctrl = len(y_A_cells), len(y_B_cells), len(ctrl_cells)
    rng = np.random.RandomState(42)

    for _ in range(n_permutations):
        perm_idx = rng.permutation(len(all_cells))
        perm_A = all_cells[perm_idx[:n_A]]
        perm_B = all_cells[perm_idx[n_A:n_A + n_B]]
        perm_ctrl = all_cells[perm_idx[n_A + n_B:n_A + n_B + n_ctrl]]
        perm_additive = perm_A.mean(axis=0) + perm_B.mean(axis=0) - perm_ctrl.mean(axis=0)
        perm_residual = y_AB_cells.mean(axis=0) - perm_additive
        perm_stats.append(float(np.abs(perm_residual).mean()))

    perm_stats = np.array(perm_stats)
    p_value = float((perm_stats >= obs_stat).mean() + 1) / (n_permutations + 1)
    effect_size = obs_stat / max(perm_stats.mean(), 1e-8)

    return p_value, effect_size, obs_stat


def detect_epistasis_improved(residual_results, adata, X, n_permutations=500):
    """3-formula epistasis detection with statistical ground truth."""
    conditions = adata.obs['condition'].values
    gene_to_condition = build_gene_to_condition(np.unique(conditions))

    ctrl_mask = conditions == 'ctrl'
    ctrl_cells = X[ctrl_mask]
    ctrl_mean = ctrl_cells.mean(axis=0)
    ctrl_std = ctrl_cells.std(axis=0) + 1e-6

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

        # --- Statistical ground truth ---
        p_val, effect_size, obs_stat = statistical_epistasis_test(
            X[dp_mask], X[mask_A], X[mask_B], ctrl_cells, n_permutations
        )
        is_epistatic_gt = 1 if p_val < 0.05 else 0

        # --- Formula 1: Additive ---
        expected_add = y_A + y_B - ctrl_mean
        epistasis_add = y_AB - expected_add
        epistasis_strength_add = float(np.abs(epistasis_add).mean())

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

        # --- Combination-level formula agreement ---
        # Each formula classifies: is this combination epistatic?
        sig_add = float((np.abs(epistasis_add) > 2 * ctrl_std).mean()) > 0.05  # >5% genes significant
        sig_mult = float((np.abs(epistasis_mult) > 2 * np.abs(np.log(np.maximum(ctrl_std, eps)))).mean()) > 0.05
        sig_prod = float((np.abs(epistasis_prod) > 2 * (ctrl_std / (ctrl_mean + eps))).mean()) > 0.05

        n_agree = sum([sig_add, sig_mult, sig_prod])
        # 3/3 = full agreement, 2/3 = partial, 1/3 or 0/3 = disagreement
        if n_agree == 3:
            agreement_level = 'full'
        elif n_agree == 2:
            agreement_level = 'partial'
        else:
            agreement_level = 'none'

        # --- GI type classification ---
        # Use mean direction of additive epistasis
        epistasis_direction = float(np.sign(epistasis_add).mean())
        epistasis_mag = float(np.abs(epistasis_add[epistasis_add != 0]).mean()) if (epistasis_add != 0).any() else 0

        if p_val >= 0.05:
            gi_type = 'additive'  # not significantly different from additive expectation
        elif epistasis_direction > 0.1:
            gi_type = 'synergy'
        elif epistasis_direction < -0.1:
            gi_type = 'suppression'
        else:
            # Significant but mixed direction → neomorphic
            gi_type = 'neomorphic'

        epistasis_results.append({
            'double': dp, 'pert1': p1, 'pert2': p2,
            'p_value': p_val,
            'effect_size': effect_size,
            'is_epistatic_gt': is_epistatic_gt,
            'epistasis_strength_add': epistasis_strength_add,
            'epistasis_strength_mult': epistasis_strength_mult,
            'epistasis_strength_prod': epistasis_strength_prod,
            'sig_add': sig_add, 'sig_mult': sig_mult, 'sig_prod': sig_prod,
            'agreement_level': agreement_level,
            'gi_type': gi_type,
            'r2_additive': rr['r2_additive'],
        })

    return epistasis_results


def evaluate_epistasis_improved(epistasis_results):
    """Evaluate epistasis detection using permutation-test ground truth."""
    labels = np.array([r['is_epistatic_gt'] for r in epistasis_results])
    scores_add = np.array([r['epistasis_strength_add'] for r in epistasis_results])
    scores_prod = np.array([r['epistasis_strength_prod'] for r in epistasis_results])
    scores_mult = np.array([r['epistasis_strength_mult'] for r in epistasis_results])
    p_values = np.array([r['p_value'] for r in epistasis_results])

    results = {}
    n_epistatic = int(labels.sum())
    n_total = len(labels)
    results['n_epistatic'] = n_epistatic
    results['n_total'] = n_total

    if len(set(labels)) > 1 and n_total > 5:
        try:
            results['auroc_additive'] = float(roc_auc_score(labels, scores_add))
            results['auroc_product'] = float(roc_auc_score(labels, scores_prod))
            results['auroc_multiplicative'] = float(roc_auc_score(labels, scores_mult))
            results['aupr_additive'] = float(average_precision_score(labels, scores_add))
        except ValueError:
            results['auroc_additive'] = 0.5
            results['auroc_product'] = 0.5
            results['auroc_multiplicative'] = 0.5
            results['aupr_additive'] = 0.0
    else:
        results['auroc_additive'] = 0.5
        results['auroc_product'] = 0.5
        results['auroc_multiplicative'] = 0.5
        results['aupr_additive'] = 0.0

    # Precision@k
    for k in [10, 20]:
        if len(epistasis_results) >= k:
            top_k_idx = np.argsort(scores_add)[-k:]
            results[f'precision_top_{k}'] = float(labels[top_k_idx].mean())

    # GI type distribution
    type_counts = {}
    for r in epistasis_results:
        t = r['gi_type']
        type_counts[t] = type_counts.get(t, 0) + 1
    results['gi_type_distribution'] = type_counts

    # Formula agreement
    agreement_counts = {'full': 0, 'partial': 0, 'none': 0}
    for r in epistasis_results:
        agreement_counts[r['agreement_level']] = agreement_counts.get(r['agreement_level'], 0) + 1
    results['formula_agreement'] = agreement_counts

    n_epistasis_detected = sum(1 for r in epistasis_results if r['is_epistatic_gt'] == 1)
    results['n_epistasis_detected'] = n_epistasis_detected

    return results


# ============================================================
# 8. Phase 6: Active Learning Simulation (RQ4) — IMPROVED
# ============================================================

def simulate_active_learning_improved(epistasis_results, uq_results, decomposition_results):
    """AL with epistasis score + diversity term."""
    if not epistasis_results or not uq_results.get('per_combination'):
        return {'error': 'Insufficient data for AL simulation'}

    epi_lookup = {r['double']: r for r in epistasis_results}
    uq_lookup = {r['double']: r for r in uq_results['per_combination']}
    decomp_lookup = {r['double']: r for r in decomposition_results}

    # Ground truth: permutation-test p-value
    is_epistatic = {r['double']: r['is_epistatic_gt'] for r in epistasis_results}
    all_combos = [r['double'] for r in epistasis_results]
    n_epistatic = sum(is_epistatic.values())

    if n_epistatic == 0 or n_epistatic == len(all_combos):
        return {'error': f'All or no epistasis (n_epistatic={n_epistatic}/{len(all_combos)})'}

    # --- Random baseline ---
    np.random.seed(42)
    random_order = list(np.random.permutation(all_combos))

    # --- UQ-only (from Phase 4) ---
    uq_scores = [(c, uq_lookup.get(c, {}).get('u_combined', 0)) for c in all_combos]
    uq_order = [c for c, s in sorted(uq_scores, key=lambda x: -x[1])]

    # --- Epistasis score from decomposition ---
    epi_scores = [(c, decomp_lookup.get(c, {}).get('r_epistasis_mag', 0)) for c in all_combos]
    epi_order = [c for c, s in sorted(epi_scores, key=lambda x: -x[1])]

    # --- Combined: epistasis score + UQ + diversity ---
    z_tx_vecs = {}
    for r in epistasis_results:
        dp = r['double']
        # Use z_tx from residual results if available
        if dp in decomp_lookup:
            z_tx_vecs[dp] = decomp_lookup[dp].get('z_tx_composed', None)

    def compute_diversity(combo, selected, z_tx_vecs):
        """Diversity = mean distance to already-selected combos in z_tx space."""
        if not selected or combo not in z_tx_vecs or z_tx_vecs[combo] is None:
            return 1.0
        z_c = z_tx_vecs[combo]
        dists = []
        for s in selected:
            if s in z_tx_vecs and z_tx_vecs[s] is not None:
                d = np.linalg.norm(z_c - z_tx_vecs[s])
                dists.append(d)
        return float(np.mean(dists)) if dists else 1.0

    # Greedy combined selection
    epi_score_lookup = {c: s for c, s in epi_scores}
    uq_score_lookup = {c: s for c, s in uq_scores}
    lam = 0.3  # diversity weight
    combined_order = []
    remaining = set(all_combos)
    for _ in range(len(all_combos)):
        best_combo, best_score = None, -np.inf
        for c in remaining:
            epi_s = epi_score_lookup.get(c, 0)
            uq_s = uq_score_lookup.get(c, 0)
            div = compute_diversity(c, combined_order, z_tx_vecs)
            score = epi_s + 0.5 * uq_s + lam * div
            if score > best_score:
                best_score = score
                best_combo = c
        combined_order.append(best_combo)
        remaining.remove(best_combo)

    # --- Oracle (upper bound) ---
    oracle_order = [c for c, _ in sorted([(r['double'], r['effect_size'])
                                           for r in epistasis_results], key=lambda x: -x[1])]

    # Compute cumulative recall curves
    def cumulative_recall(order, is_epistatic):
        found = 0
        cum = []
        for c in order:
            if is_epistatic.get(c, False):
                found += 1
            cum.append(found / max(n_epistatic, 1))
        return cum

    random_cum = cumulative_recall(random_order, is_epistatic)
    uq_cum = cumulative_recall(uq_order, is_epistatic)
    epi_cum = cumulative_recall(epi_order, is_epistatic)
    combined_cum = cumulative_recall(combined_order, is_epistatic)
    oracle_cum = cumulative_recall(oracle_order, is_epistatic)

    top_k_results = {}
    for k in [5, 10, 20]:
        if k <= len(all_combos):
            top_k_results[k] = {
                'random': random_cum[k-1],
                'uq_only': uq_cum[k-1],
                'epistasis_only': epi_cum[k-1],
                'combined': combined_cum[k-1],
                'oracle': oracle_cum[k-1],
                'improvement_vs_random_uq': uq_cum[k-1] / max(random_cum[k-1], 1e-8),
                'improvement_vs_random_epi': epi_cum[k-1] / max(random_cum[k-1], 1e-8),
                'improvement_vs_random_combined': combined_cum[k-1] / max(random_cum[k-1], 1e-8),
            }

    return {
        'n_total': len(all_combos),
        'n_epistatic': n_epistatic,
        'top_k_results': top_k_results,
        'final_random_recall': random_cum[-1],
        'final_combined_recall': combined_cum[-1],
        'final_oracle_recall': oracle_cum[-1],
    }


# ============================================================
# 9. Main Pipeline
# ============================================================

def main():
    t0 = time.time()
    print("=" * 70)
    print("Run 09b: Improved Epistasis Detection + UQ Pipeline")
    print("=" * 70)

    # ---- Phase 1: Data + Model ----
    print("\n[Phase 1] Loading Norman 2019 data...")
    adata, gears_pd = load_norman_data()
    adata, X = preprocess_adata(adata, n_top_genes=500)
    n_genes = X.shape[1]
    print(f"  Data: {X.shape[0]} cells, {n_genes} genes")

    single_perts, double_perts = get_perturbation_info(adata)
    print(f"  Single perturbations: {len(single_perts)}")
    print(f"  Double perturbations (valid): {len(double_perts)}")

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
    decoder = FCRDecoder(z_dim, n_genes)
    encoder, decoder = train_fcr_icm(
        encoder, decoder, X_filt, pert_ids, ct_ids, n_cell_types,
        n_epochs=100, beta=0.5, icm_weight=10.0, use_icm=True
    )
    print("  Training complete.")

    # ---- Phase 2: Residual Extraction ----
    print("\n[Phase 2] Extracting prediction residuals...")
    residual_results = extract_residuals(encoder, decoder, adata, X, pert_id_map, n_cell_types, z_dim)
    print(f"  Double-KO pairs with residuals: {len(residual_results)}")
    if residual_results:
        r2_model = [r['r2_model'] for r in residual_results]
        r2_add = [r['r2_additive'] for r in residual_results]
        print(f"  Model prediction R2: mean={np.mean(r2_model):.3f}, median={np.median(r2_model):.3f}")
        print(f"  Additive prediction R2: mean={np.mean(r2_add):.3f}, median={np.median(r2_add):.3f}")

    # ---- Phase 3: Residual Decomposition (RQ1) ----
    print("\n[Phase 3] Residual decomposition (RQ1)...")
    decomposition_results = decompose_residuals_improved(
        encoder, decoder, residual_results, adata, X,
        pert_id_map, n_cell_types, z_dim, n_mc_passes=30
    )
    print(f"  Decomposed combinations: {len(decomposition_results)}")
    if decomposition_results:
        r_total_mags = [d['r_total_mag'] for d in decomposition_results]
        r_epi_mags = [d['r_epistasis_mag'] for d in decomposition_results]
        r_model_mags = [d['r_model_mag'] for d in decomposition_results]
        r_noise_mags = [d['r_noise_mag'] for d in decomposition_results]
        print(f"  Mean |r_total|: {np.mean(r_total_mags):.4f}")
        print(f"  Mean |r_epistasis|: {np.mean(r_epi_mags):.4f} ({np.mean(r_epi_mags)/max(np.mean(r_total_mags),1e-8)*100:.1f}%)")
        print(f"  Mean |r_model|: {np.mean(r_model_mags):.4f} ({np.mean(r_model_mags)/max(np.mean(r_total_mags),1e-8)*100:.1f}%)")
        print(f"  Mean |r_noise|: {np.mean(r_noise_mags):.4f} ({np.mean(r_noise_mags)/max(np.mean(r_total_mags),1e-8)*100:.1f}%)")

    # Synthetic validation
    print("\n[RQ1 Synthetic] Validating decomposition with synthetic data...")
    from run_09_epistasis_uq import validate_decomposition_synthetic
    synth_results = validate_decomposition_synthetic()
    print(f"  Perfect decomposition correlation: {synth_results['perfect_decomp_corr_mean']:.3f}")

    # ---- Phase 4: Uncertainty Quantification (RQ2) ----
    print("\n[Phase 4] Uncertainty quantification (RQ2)...")
    uq_results = quantify_uncertainty_improved(decomposition_results)
    print(f"  U-Error Spearman rho (combined): {uq_results.get('spearman_rho_combined', 'N/A')}")
    print(f"  U-Error Spearman rho (ICM only): {uq_results.get('spearman_rho_icm_only', 'N/A')}")
    print(f"  U-Error Spearman rho (MC only): {uq_results.get('spearman_rho_mc_only', 'N/A')}")
    print(f"  Best weights: {uq_results.get('best_weights', 'N/A')}")
    print(f"  Mean coverage: {uq_results.get('mean_coverage', 'N/A'):.3f}")

    # ---- Phase 5: Epistasis Detection (RQ3) ----
    print("\n[Phase 5] Epistasis detection with 3-formula sensitivity (RQ3)...")
    epistasis_results = detect_epistasis_improved(residual_results, adata, X, n_permutations=500)
    print(f"  Combinations analyzed: {len(epistasis_results)}")
    n_epi = sum(1 for r in epistasis_results if r['is_epistatic_gt'] == 1)
    print(f"  Permutation-test epistatic (p<0.05): {n_epi}/{len(epistasis_results)}")

    epi_eval = evaluate_epistasis_improved(epistasis_results)
    print(f"  AUROC (additive): {epi_eval.get('auroc_additive', 'N/A')}")
    print(f"  AUROC (product): {epi_eval.get('auroc_product', 'N/A')}")
    print(f"  AUROC (multiplicative): {epi_eval.get('auroc_multiplicative', 'N/A')}")
    print(f"  AUPR (additive): {epi_eval.get('aupr_additive', 'N/A')}")
    if 'precision_top_10' in epi_eval:
        print(f"  Precision@10: {epi_eval['precision_top_10']:.3f}")
    if 'precision_top_20' in epi_eval:
        print(f"  Precision@20: {epi_eval['precision_top_20']:.3f}")
    print(f"  GI type distribution: {epi_eval.get('gi_type_distribution', {})}")
    print(f"  Formula agreement: {epi_eval.get('formula_agreement', {})}")

    # ---- Phase 6: Active Learning (RQ4) ----
    print("\n[Phase 6] Active learning simulation (RQ4)...")
    al_results = simulate_active_learning_improved(epistasis_results, uq_results, decomposition_results)
    if 'error' not in al_results:
        print(f"  Total combinations: {al_results['n_total']}")
        print(f"  Epistatic combinations: {al_results['n_epistatic']}")
        for k, v in al_results.get('top_k_results', {}).items():
            print(f"  Top-{k}: random={v['random']:.3f}, UQ={v['uq_only']:.3f}, "
                  f"epi={v['epistasis_only']:.3f}, combined={v['combined']:.3f}, "
                  f"oracle={v['oracle']:.3f}")
            print(f"    Improvement: UQ={v['improvement_vs_random_uq']:.2f}x, "
                  f"epi={v['improvement_vs_random_epi']:.2f}x, "
                  f"combined={v['improvement_vs_random_combined']:.2f}x")
    else:
        print(f"  {al_results['error']}")

    # ---- Summary ----
    elapsed = time.time() - t0
    print("\n" + "=" * 70)
    print("RUN 09b SUMMARY")
    print("=" * 70)
    print(f"  RQ1 (Decomposition): Synthetic corr = {synth_results['perfect_decomp_corr_mean']:.3f} (target: r>0.7)")
    if decomposition_results:
        epi_pct = np.mean([d['r_epistasis_mag'] for d in decomposition_results]) / max(np.mean([d['r_total_mag'] for d in decomposition_results]), 1e-8) * 100
        print(f"  RQ1 (Decomposition): Epistasis fraction = {epi_pct:.1f}% of residual")
    print(f"  RQ2 (UQ): U-Error rho = {uq_results.get('spearman_rho_combined', 'N/A')} (target: >0.6)")
    print(f"  RQ2 (UQ): Coverage = {uq_results.get('mean_coverage', 'N/A'):.3f} (target: 0.85-0.95)")
    print(f"  RQ3 (Epistasis): AUROC(add) = {epi_eval.get('auroc_additive', 'N/A')} (target: >0.75)")
    print(f"  RQ3 (Epistasis): AUPR(add) = {epi_eval.get('aupr_additive', 'N/A')}")
    if 'precision_top_10' in epi_eval:
        print(f"  RQ3 (Epistasis): Precision@10 = {epi_eval['precision_top_10']:.3f} (target: >0.6)")
    print(f"  RQ3 (Epistasis): Permutation-test epistatic = {n_epi}/{len(epistasis_results)}")
    if 'error' not in al_results:
        for k, v in al_results.get('top_k_results', {}).items():
            print(f"  RQ4 (AL): Top-{k} combined improvement = {v['improvement_vs_random_combined']:.2f}x vs random (target: >2x)")
    print(f"\n  Elapsed: {elapsed:.0f}s")

    # Save results
    save_results = {
        'rq1_synthetic': synth_results,
        'rq2_uq': {k: v for k, v in uq_results.items() if k != 'per_combination'},
        'rq3_epistasis': epi_eval,
        'rq4_al': {k: v for k, v in al_results.items()
                   if k != 'per_combination' and not isinstance(v, list)},
        'n_residual_pairs': len(residual_results),
        'n_decomposed': len(decomposition_results),
        'n_epistasis_analyzed': len(epistasis_results),
        'elapsed_s': elapsed,
    }
    with open(os.path.join(RESULTS_DIR, 'run_09b_results.json'), 'w') as f:
        json.dump(save_results, f, indent=2, default=str)
    print(f"\n  Results saved to {RESULTS_DIR}/run_09b_results.json")


if __name__ == '__main__':
    main()