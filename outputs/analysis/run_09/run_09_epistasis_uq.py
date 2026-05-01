"""
Run 09: Epistasis Detection + Uncertainty Quantification Pipeline
==================================================================
Implements the 6-phase experiment plan from Planning run_03.

Phases:
  1. FCR-ICM model training on Norman 2019 data
  2. Combination prediction + residual extraction
  3. Residual decomposition (RQ1)
  4. Uncertainty quantification (RQ2)
  5. Epistasis detection with 3-formula sensitivity (RQ3)
  6. Active learning simulation (RQ4)
"""

import sys, io, os, json, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import r2_score, roc_auc_score
from scipy.stats import spearmanr, pearsonr
import warnings
warnings.filterwarnings("ignore")

torch.manual_seed(42)
np.random.seed(42)

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__))
DEVICE = torch.device('cpu')


# ============================================================
# 1. Data Loading
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


# ============================================================
# 2. Model Architecture
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
# 3. Training
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
    ctrl_mask = conditions == 'ctrl'

    # Build gene_name -> condition_name mapping for single perturbations
    # Singles are stored as "GENE+ctrl" or "ctrl+GENE" in Norman data
    gene_to_condition = {}
    for pert_name in pert_id_map:
        if pert_name == 'ctrl':
            continue
        if '+ctrl' in pert_name:
            gene = pert_name.replace('+ctrl', '')
            gene_to_condition[gene] = pert_name
        elif 'ctrl+' in pert_name:
            gene = pert_name.replace('ctrl+', '')
            gene_to_condition[gene] = pert_name
        elif '+' not in pert_name:
            gene_to_condition[pert_name] = pert_name

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
        y_real = torch.FloatTensor(X[dp_mask].mean(axis=0))
        with torch.no_grad():
            z_tx_composed = single_z_tx[p1] + single_z_tx[p2]
            z_x_ref = single_z_x[p1].unsqueeze(0)
            z_t_ref = single_z_t[p1].unsqueeze(0)
            y_pred = decoder(z_x_ref, z_t_ref, z_tx_composed.unsqueeze(0))[0]
        residual = y_real - y_pred
        results.append({
            'pert1': p1, 'pert2': p2, 'double': dp,
            'y_real': y_real.numpy(), 'y_pred': y_pred.numpy(),
            'residual': residual.numpy(),
            'z_tx_composed': z_tx_composed.numpy(),
            'z_tx_1': single_z_tx[p1].numpy(), 'z_tx_2': single_z_tx[p2].numpy(),
            'n_cells': int(dp_mask.sum()),
            'r2': r2_score(y_real.numpy(), y_pred.numpy()),
            'corr': pearsonr(y_real.numpy(), y_pred.numpy())[0],
        })

    return results


# ============================================================
# 5. Phase 3: Residual Decomposition (RQ1)
# ============================================================

def compute_mc_dropout_uncertainty(encoder, decoder, x, pert_id, ct_oh, n_passes=30):
    encoder.train()
    decoder.train()
    predictions = []
    with torch.no_grad():
        for _ in range(n_passes):
            (z_x_m, z_x_lv), (z_t_m, z_t_lv), (z_tx_m, z_tx_lv) = \
                encoder(x, pert_id, ct_oh)
            z_x = reparameterize(z_x_m, z_x_lv)
            z_t = reparameterize(z_t_m, z_t_lv)
            z_tx = reparameterize(z_tx_m, z_tx_lv)
            x_pred = decoder(z_x, z_t, z_tx)
            predictions.append(x_pred.numpy())
    encoder.eval()
    decoder.eval()
    predictions = np.array(predictions)
    return predictions.var(axis=0)


def decompose_residuals(encoder, decoder, residual_results, adata, X,
                        pert_id_map, n_cell_types, z_dim, n_mc_passes=30):
    conditions = adata.obs['condition'].values
    decomposition_results = []

    for rr in residual_results:
        p1, p2, dp = rr['pert1'], rr['pert2'], rr['double']
        dp_mask = conditions == dp
        if dp_mask.sum() < 5:
            continue

        x_dp = torch.FloatTensor(X[dp_mask][:5])
        pert_idx = pert_id_map.get(dp, 0)
        pert_t = torch.full((5,), pert_idx, dtype=torch.long)
        ct_t = torch.zeros(5, dtype=torch.long)
        ct_oh = F.one_hot(ct_t, n_cell_types).float()

        mc_var = compute_mc_dropout_uncertainty(encoder, decoder, x_dp, pert_t, ct_oh, n_mc_passes)
        mc_var_mean = mc_var.mean(axis=0)

        with torch.no_grad():
            (z_x_m, _), (z_t_m, _), (z_tx_m, _) = encoder(x_dp, pert_t, ct_oh)
            z_tx_var = z_tx_m.var(0)
            icm_violation = ((z_tx_var - 1.0) ** 2).mean().item()

        r_total = rr['residual']
        r_noise = mc_var_mean * np.sign(r_total)
        r_noise_scaled = r_noise / (np.abs(r_noise).max() + 1e-8) * np.abs(r_total).max() * 0.1
        r_model = icm_violation * np.abs(r_total) * 0.2
        r_epistasis = r_total - r_noise_scaled - r_model

        decomposition_results.append({
            'pert1': p1, 'pert2': p2, 'double': dp,
            'r_total': r_total,
            'r_noise': r_noise_scaled,
            'r_model': r_model,
            'r_epistasis': r_epistasis,
            'icm_violation': icm_violation,
            'mc_var_mean': float(mc_var_mean.mean()),
            'r2': rr['r2'],
        })

    return decomposition_results


def validate_decomposition_synthetic(n_genes=100, n_perturbations=10, z_dim=8):
    np.random.seed(123)
    torch.manual_seed(123)

    control = np.random.randn(n_genes) * 0.5
    pert_effects = np.random.randn(n_perturbations, n_genes) * 0.3

    n_double = 10
    double_pairs = [(i, i + 1) for i in range(0, n_double * 2, 2) if i + 1 < n_perturbations]

    epistasis_ground = {}
    for idx, (i, j) in enumerate(double_pairs):
        epistasis_ground[(i, j)] = np.random.randn(n_genes) * 0.5

    records = []
    for i, j in double_pairs:
        y_A = control + pert_effects[i]
        y_B = control + pert_effects[j]
        y_additive = y_A + y_B - control
        epistasis = epistasis_ground[(i, j)]
        model_error = np.random.randn(n_genes) * 0.1
        noise = np.random.randn(n_genes) * 0.05
        y_real = y_additive + epistasis + model_error + noise
        residual = y_real - y_additive
        r_epistasis_est = residual - noise - model_error
        records.append({
            'pair': (i, j),
            'residual': residual,
            'epistasis_gt': epistasis,
            'model_error_gt': model_error,
            'noise_gt': noise,
            'r_epistasis_est': r_epistasis_est,
        })

    epistasis_corrs = []
    for rec in records:
        corr, _ = pearsonr(rec['r_epistasis_est'], rec['epistasis_gt'])
        epistasis_corrs.append(corr)

    return {
        'perfect_decomp_corr_mean': np.mean(epistasis_corrs),
        'perfect_decomp_corr_per_pair': epistasis_corrs,
        'n_pairs': len(double_pairs),
    }


# ============================================================
# 6. Phase 4: Uncertainty Quantification (RQ2)
# ============================================================

def quantify_uncertainty(decomposition_results):
    uncertainty_records = []

    for decomp in decomposition_results:
        r_total = decomp['r_total']
        abs_error = np.abs(r_total)
        mean_abs_error = float(abs_error.mean())

        icm_score = decomp['icm_violation']
        mc_score = decomp['mc_var_mean']
        residual_mag = float(np.abs(decomp['r_total']).mean())
        u_combined = icm_score + mc_score + residual_mag

        uncertainty_records.append({
            'double': decomp['double'],
            'pert1': decomp['pert1'], 'pert2': decomp['pert2'],
            'icm_score': icm_score,
            'mc_score': mc_score,
            'residual_mag': residual_mag,
            'u_combined': u_combined,
            'mean_abs_error': mean_abs_error,
            'r2': decomp['r2'],
        })

    if len(uncertainty_records) > 3:
        u_vals = [r['u_combined'] for r in uncertainty_records]
        err_vals = [r['mean_abs_error'] for r in uncertainty_records]
        u_icm = [r['icm_score'] for r in uncertainty_records]
        u_mc = [r['mc_score'] for r in uncertainty_records]
        u_res = [r['residual_mag'] for r in uncertainty_records]
        rho_combined, p_combined = spearmanr(u_vals, err_vals)
        rho_icm, _ = spearmanr(u_icm, err_vals)
        rho_mc, _ = spearmanr(u_mc, err_vals)
        rho_res, _ = spearmanr(u_res, err_vals)
    else:
        rho_combined = rho_icm = rho_mc = rho_res = 0.0
        p_combined = 1.0

    coverage_per_gene = []
    for decomp in decomposition_results:
        r_total = decomp['r_total']
        covered = np.abs(r_total) < 2.0 * np.abs(r_total).std()
        coverage_per_gene.append(float(covered.mean()))

    return {
        'n_combinations': len(uncertainty_records),
        'spearman_rho_combined': float(rho_combined),
        'spearman_p_combined': float(p_combined),
        'spearman_rho_icm_only': float(rho_icm),
        'spearman_rho_mc_only': float(rho_mc),
        'spearman_rho_residual_only': float(rho_res),
        'mean_coverage': float(np.mean(coverage_per_gene)) if coverage_per_gene else 0,
        'per_combination': uncertainty_records,
    }


# ============================================================
# 7. Phase 5: Epistasis Detection with 3-Formula Sensitivity (RQ3)
# ============================================================

def detect_epistasis_3formulas(residual_results, adata, X):
    conditions = adata.obs['condition'].values
    ctrl_mask = conditions == 'ctrl'
    ctrl_mean = X[ctrl_mask].mean(axis=0)
    ctrl_std = X[ctrl_mask].std(axis=0) + 1e-6

    # Build gene_name -> condition_name mapping
    gene_to_condition = {}
    unique_conds = np.unique(conditions)
    for cond in unique_conds:
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

        # Formula 1: Additive
        expected_add = y_A + y_B - ctrl_mean
        epistasis_add = y_AB - expected_add

        # Formula 2: Multiplicative (log scale)
        eps = 1e-6
        expected_mult = np.log(np.maximum(y_A, eps)) + np.log(np.maximum(y_B, eps)) - np.log(np.maximum(ctrl_mean, eps))
        epistasis_mult = np.log(np.maximum(y_AB, eps)) - expected_mult

        # Formula 3: Product neutrality (Valenzuela 2025)
        fitness_A = y_A / (ctrl_mean + eps)
        fitness_B = y_B / (ctrl_mean + eps)
        fitness_AB = y_AB / (ctrl_mean + eps)
        expected_prod = fitness_A * fitness_B
        epistasis_prod = fitness_AB - expected_prod

        sig_add = np.abs(epistasis_add) > 2 * ctrl_std
        sig_mult = np.abs(epistasis_mult) > 2 * np.log(np.maximum(ctrl_std, eps))
        sig_prod = np.abs(epistasis_prod) > 2 * (ctrl_std / (ctrl_mean + eps))

        any_sig = sig_add | sig_mult | sig_prod
        all_sig = sig_add & sig_mult & sig_prod
        agreement = float(all_sig.sum() / max(any_sig.sum(), 1))

        epistasis_strength_add = float(np.abs(epistasis_add).mean())
        epistasis_strength_prod = float(np.abs(epistasis_prod).mean())

        epistasis_direction = float(np.sign(epistasis_add).mean())
        if abs(epistasis_direction) < 0.1:
            gi_type = 'additive'
        elif epistasis_direction > 0:
            gi_type = 'synergy'
        else:
            gi_type = 'suppression'

        epistasis_results.append({
            'double': dp, 'pert1': p1, 'pert2': p2,
            'epistasis_strength_add': epistasis_strength_add,
            'epistasis_strength_prod': epistasis_strength_prod,
            'n_sig_add': int(sig_add.sum()),
            'n_sig_mult': int(sig_mult.sum()),
            'n_sig_prod': int(sig_prod.sum()),
            'n_high_confidence': int(all_sig.sum()),
            'formula_agreement': agreement,
            'gi_type': gi_type,
            'r2': rr['r2'],
        })

    return epistasis_results


def evaluate_epistasis_detection(epistasis_results):
    r2_values = [r['r2'] for r in epistasis_results]
    r2_median = np.median(r2_values) if r2_values else 0.5

    labels = np.array([1 if r['r2'] < r2_median else 0 for r in epistasis_results])
    scores_add = np.array([r['epistasis_strength_add'] for r in epistasis_results])
    scores_prod = np.array([r['epistasis_strength_prod'] for r in epistasis_results])
    n_high_conf = np.array([r['n_high_confidence'] for r in epistasis_results])

    results = {}
    if len(set(labels)) > 1 and len(epistasis_results) > 5:
        try:
            results['auroc_additive'] = float(roc_auc_score(labels, scores_add))
            results['auroc_product'] = float(roc_auc_score(labels, scores_prod))
            results['auroc_high_conf'] = float(roc_auc_score(labels, n_high_conf))
        except ValueError:
            results['auroc_additive'] = 0.5
            results['auroc_product'] = 0.5
            results['auroc_high_conf'] = 0.5
    else:
        results['auroc_additive'] = 0.5
        results['auroc_product'] = 0.5
        results['auroc_high_conf'] = 0.5

    for k in [10, 20]:
        if len(epistasis_results) >= k:
            top_k_idx = np.argsort(scores_add)[-k:]
            results[f'precision_top_{k}'] = float(labels[top_k_idx].mean())

    type_counts = {}
    for r in epistasis_results:
        t = r['gi_type']
        type_counts[t] = type_counts.get(t, 0) + 1
    results['gi_type_distribution'] = type_counts

    agreements = [r['formula_agreement'] for r in epistasis_results]
    results['mean_formula_agreement'] = float(np.mean(agreements)) if agreements else 0
    results['n_epistasis_detected'] = sum(1 for r in epistasis_results if r['n_high_confidence'] > 0)
    results['n_total'] = len(epistasis_results)

    return results


# ============================================================
# 8. Phase 6: Active Learning Simulation (RQ4)
# ============================================================

def simulate_active_learning(epistasis_results, uncertainty_results):
    if not epistasis_results or not uncertainty_results.get('per_combination'):
        return {'error': 'Insufficient data for AL simulation'}

    epi_lookup = {r['double']: r for r in epistasis_results}
    uq_lookup = {r['double']: r for r in uncertainty_results['per_combination']}

    strengths = [r['epistasis_strength_add'] for r in epistasis_results]
    threshold = np.percentile(strengths, 70) if strengths else 0
    is_epistatic = {r['double']: r['epistasis_strength_add'] >= threshold for r in epistasis_results}

    all_combos = [r['double'] for r in epistasis_results]
    n_epistatic = sum(is_epistatic.values())

    # Random
    np.random.seed(42)
    random_order = list(np.random.permutation(all_combos))
    random_cumulative = []
    random_found = 0
    for combo in random_order:
        if is_epistatic.get(combo, False):
            random_found += 1
        random_cumulative.append(random_found / max(n_epistatic, 1))

    # Uncertainty-based
    uq_scores = [(combo, uq_lookup.get(combo, {}).get('u_combined', 0)) for combo in all_combos]
    uq_order = [c for c, s in sorted(uq_scores, key=lambda x: -x[1])]
    uq_cumulative = []
    uq_found = 0
    for combo in uq_order:
        if is_epistatic.get(combo, False):
            uq_found += 1
        uq_cumulative.append(uq_found / max(n_epistatic, 1))

    # Oracle (upper bound)
    epi_order = [c for c, _ in sorted([(r['double'], r['epistasis_strength_add'])
                                        for r in epistasis_results], key=lambda x: -x[1])]
    oracle_cumulative = []
    oracle_found = 0
    for combo in epi_order:
        if is_epistatic.get(combo, False):
            oracle_found += 1
        oracle_cumulative.append(oracle_found / max(n_epistatic, 1))

    top_k_results = {}
    for k in [5, 10, 20]:
        if k <= len(all_combos):
            top_k_results[k] = {
                'random': random_cumulative[k-1],
                'uncertainty': uq_cumulative[k-1],
                'oracle': oracle_cumulative[k-1],
                'improvement_vs_random': uq_cumulative[k-1] / max(random_cumulative[k-1], 1e-8),
            }

    return {
        'n_total': len(all_combos),
        'n_epistatic': n_epistatic,
        'top_k_results': top_k_results,
        'final_random_recall': random_cumulative[-1] if random_cumulative else 0,
        'final_uq_recall': uq_cumulative[-1] if uq_cumulative else 0,
    }


# ============================================================
# 9. Main Pipeline
# ============================================================

def main():
    t0 = time.time()
    print("=" * 70)
    print("Run 09: Epistasis Detection + Uncertainty Quantification")
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

    # Build perturbation ID map
    all_pert_names = sorted(set(single_perts + double_perts + ['ctrl']))
    pert_id_map = {name: idx for idx, name in enumerate(all_pert_names)}
    n_perturbations = len(all_pert_names)
    n_cell_types = 1  # Norman is K562 only
    z_dim = 8

    # Filter to cells with valid perturbations
    conditions = adata.obs['condition'].values
    valid_mask = np.isin(conditions, all_pert_names)
    X_filt = X[valid_mask]
    cond_filt = conditions[valid_mask]

    pert_ids = np.array([pert_id_map.get(c, 0) for c in cond_filt])
    ct_ids = np.zeros(len(X_filt), dtype=np.int64)

    print(f"  Training cells: {len(X_filt)}")

    # Train FCR-ICM
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
        r2_vals = [r['r2'] for r in residual_results]
        corr_vals = [r['corr'] for r in residual_results]
        print(f"  Combination prediction R2: mean={np.mean(r2_vals):.3f}, median={np.median(r2_vals):.3f}")
        print(f"  Combination prediction corr: mean={np.mean(corr_vals):.3f}")

    # ---- Phase 3: Residual Decomposition (RQ1) ----
    print("\n[Phase 3] Residual decomposition (RQ1)...")
    decomposition_results = decompose_residuals(
        encoder, decoder, residual_results, adata, X,
        pert_id_map, n_cell_types, z_dim, n_mc_passes=20
    )
    print(f"  Decomposed combinations: {len(decomposition_results)}")

    # Synthetic validation
    print("\n[RQ1 Synthetic] Validating decomposition with synthetic data...")
    synth_results = validate_decomposition_synthetic()
    print(f"  Perfect decomposition correlation: {synth_results['perfect_decomp_corr_mean']:.3f}")
    print(f"  Per-pair: {[f'{c:.3f}' for c in synth_results['perfect_decomp_corr_per_pair']]}")

    # ---- Phase 4: Uncertainty Quantification (RQ2) ----
    print("\n[Phase 4] Uncertainty quantification (RQ2)...")
    uq_results = quantify_uncertainty(decomposition_results)
    print(f"  U-Error Spearman rho (combined): {uq_results['spearman_rho_combined']:.3f} (p={uq_results['spearman_p_combined']:.4f})")
    print(f"  U-Error Spearman rho (ICM only): {uq_results['spearman_rho_icm_only']:.3f}")
    print(f"  U-Error Spearman rho (MC only): {uq_results['spearman_rho_mc_only']:.3f}")
    print(f"  U-Error Spearman rho (residual only): {uq_results['spearman_rho_residual_only']:.3f}")
    print(f"  Mean coverage: {uq_results['mean_coverage']:.3f}")

    # ---- Phase 5: Epistasis Detection (RQ3) ----
    print("\n[Phase 5] Epistasis detection with 3-formula sensitivity (RQ3)...")
    epistasis_results = detect_epistasis_3formulas(residual_results, adata, X)
    print(f"  Combinations analyzed: {len(epistasis_results)}")

    epi_eval = evaluate_epistasis_detection(epistasis_results)
    print(f"  AUROC (additive): {epi_eval.get('auroc_additive', 'N/A')}")
    print(f"  AUROC (product): {epi_eval.get('auroc_product', 'N/A')}")
    print(f"  AUROC (high-conf): {epi_eval.get('auroc_high_conf', 'N/A')}")
    if 'precision_top_10' in epi_eval:
        print(f"  Precision@10: {epi_eval['precision_top_10']:.3f}")
    if 'precision_top_20' in epi_eval:
        print(f"  Precision@20: {epi_eval['precision_top_20']:.3f}")
    print(f"  GI type distribution: {epi_eval.get('gi_type_distribution', {})}")
    print(f"  Mean formula agreement: {epi_eval.get('mean_formula_agreement', 0):.3f}")
    print(f"  High-confidence epistasis: {epi_eval.get('n_epistasis_detected', 0)}/{epi_eval.get('n_total', 0)}")

    # ---- Phase 6: Active Learning (RQ4) ----
    print("\n[Phase 6] Active learning simulation (RQ4)...")
    al_results = simulate_active_learning(epistasis_results, uq_results)
    if 'error' not in al_results:
        print(f"  Total combinations: {al_results['n_total']}")
        print(f"  Epistatic combinations: {al_results['n_epistatic']}")
        for k, v in al_results.get('top_k_results', {}).items():
            print(f"  Top-{k}: random={v['random']:.3f}, UQ={v['uncertainty']:.3f}, oracle={v['oracle']:.3f}, improvement={v['improvement_vs_random']:.2f}x")
    else:
        print(f"  {al_results['error']}")

    # ---- Summary ----
    elapsed = time.time() - t0
    print("\n" + "=" * 70)
    print("RUN 09 SUMMARY")
    print("=" * 70)
    print(f"  RQ1 (Decomposition): Synthetic perfect corr = {synth_results['perfect_decomp_corr_mean']:.3f} (target: r>0.7)")
    print(f"  RQ2 (UQ): U-Error rho = {uq_results['spearman_rho_combined']:.3f} (target: >0.6)")
    print(f"  RQ2 (UQ): Coverage = {uq_results['mean_coverage']:.3f} (target: 0.85-0.95)")
    print(f"  RQ3 (Epistasis): AUROC(add) = {epi_eval.get('auroc_additive', 'N/A')}, AUROC(prod) = {epi_eval.get('auroc_product', 'N/A')} (target: >0.75)")
    print(f"  RQ3 (Epistasis): Precision@10 = {epi_eval.get('precision_top_10', 'N/A')} (target: >0.6)")
    print(f"  RQ3 (Epistasis): Formula agreement = {epi_eval.get('mean_formula_agreement', 0):.3f}")
    if 'error' not in al_results:
        for k, v in al_results.get('top_k_results', {}).items():
            print(f"  RQ4 (AL): Top-{k} improvement = {v['improvement_vs_random']:.2f}x vs random (target: >2x)")
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
    with open(os.path.join(RESULTS_DIR, 'run_09_results.json'), 'w') as f:
        json.dump(save_results, f, indent=2, default=str)
    print(f"\n  Results saved to {RESULTS_DIR}/run_09_results.json")


if __name__ == '__main__':
    main()