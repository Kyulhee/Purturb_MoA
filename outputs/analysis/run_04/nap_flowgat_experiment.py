"""
NAP FlowGAT-style Experiment: Condition-specific Graph Variation (C2+C6)
=========================================================================
Core insight: The original experiment failed because all knockout samples
share the SAME graph topology (C2=F). If we make the graph VARY across
samples by incorporating FBA-derived edge weights, then:
  - C2=T: graphs differ per knockout condition (Mass Flow Graphs)
  - C6=T: edge weights carry flux/flow information

This is the FlowGAT approach (Hasibi et al., 2024, NPJ Syst Biol Appl):
  - For each knockout condition, run FBA to get flux distribution
  - Construct a Mass Flow Graph (MFG) where edge weights = flux values
  - Feed condition-specific MFG to GNN → different graph per sample

NAP Prediction: C2=T + C6=T → 2/6 → GNN value expected

EXPERIMENT DESIGN (honest test):
  A. XGBoost-only baseline: knockout mask → growth (no graph)
  B. GNN+XGBoost (fixed graph): same graph for all samples (C2=F, C6=F)
  C. GNN+XGBoost (MFG): different graph per sample (C2=T, C6=T)

  If NAP is correct: C > A, but B ≈ A or B < A

CRITICAL NOTE: There is a circular logic problem — we need FBA to create MFG,
but the surrogate's purpose is to REPLACE FBA. This experiment tests WHETHER
MFG-based GNN provides value, not whether it's practically deployable.
The practical solution would be: coarse FBA/pFBA for MFG, fine surrogate for growth.

Micro-validation: textbook model, 200 samples, 3-fold CV
"""

import time
import random
import warnings
from typing import Dict, List, Optional, Tuple

import cobra
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import HeteroData
from torch_geometric.nn import HGTConv, GATConv, GCNConv
import xgboost as xgb
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.feature_selection import mutual_info_regression

warnings.filterwarnings("ignore")


# ============================================================
# 1. Model & Data (reuse from nap_transfer_experiment)
# ============================================================

def load_textbook():
    model = cobra.io.load_model("textbook")
    return model


def run_fba_knockout(model, gene_ids):
    try:
        model_cp = model.copy()
        with model_cp:
            for gid in gene_ids:
                if gid in [g.id for g in model_cp.genes]:
                    model_cp.genes.get_by_id(gid).knock_out()
            sol = model_cp.optimize()
            if sol.status == "optimal":
                return sol.objective_value, sol.fluxes
        return 0.0, None
    except:
        return 0.0, None


def generate_random_combos(model, n, min_k=1, max_k=5, seed=42):
    gene_ids = [g.id for g in model.genes]
    rng = random.Random(seed)
    combos = []
    for _ in range(n):
        k = rng.randint(min_k, min(max_k, len(gene_ids)))
        combo = rng.sample(gene_ids, k)
        combos.append(combo)
    return combos


def build_knockout_mask(model, combos):
    gene_ids = sorted([g.id for g in model.genes])
    gene_idx = {gid: i for i, gid in enumerate(gene_ids)}
    n_gene = len(gene_ids)
    mask = np.zeros((len(combos), n_gene), dtype=np.float32)
    for i, combo in enumerate(combos):
        for gid in combo:
            if gid in gene_idx:
                mask[i, gene_idx[gid]] = 1.0
    return mask


# ============================================================
# 2. Mass Flow Graph Construction (FlowGAT-style)
# ============================================================

def build_mfg_for_knockout(model, gene_ids, knockout_genes):
    """
    Build Mass Flow Graph for a specific knockout condition.

    After knocking out genes, run FBA to get flux distribution.
    Then construct a HeteroData graph where edge weights = flux values.

    This makes the graph DIFFERENT for each knockout condition (C2=T),
    and edge weights carry physical information (C6=T).
    """
    model_cp = model.copy()

    # Apply knockouts
    with model_cp:
        for gid in knockout_genes:
            if gid in [g.id for g in model_cp.genes]:
                model_cp.genes.get_by_id(gid).knock_out()

        sol = model_cp.optimize()
        if sol.status != "optimal":
            return None, 0.0

        growth = sol.objective_value
        fluxes = sol.fluxes

    # Build graph with flux-weighted edges
    data = HeteroData()

    met_ids = sorted([m.id for m in model.metabolites])
    rxn_ids = sorted([r.id for r in model.reactions])
    gene_ids_sorted = sorted([g.id for g in model.genes])

    met_idx = {mid: i for i, mid in enumerate(met_ids)}
    rxn_idx = {rid: i for i, rid in enumerate(rxn_ids)}
    gene_idx = {gid: i for i, gid in enumerate(gene_ids_sorted)}

    n_met, n_rxn, n_gene = len(met_ids), len(rxn_ids), len(gene_ids_sorted)

    # Stoichiometry edges with FLUX WEIGHTS
    consumes_src, consumes_dst, consumes_flux = [], [], []
    produces_src, produces_dst, produces_flux = [], [], []

    for rxn in model.reactions:
        ri = rxn_idx[rxn.id]
        flux_val = abs(fluxes[rxn.id]) if rxn.id in fluxes.index else 0.0

        for met, coeff in rxn.metabolites.items():
            mi = met_idx[met.id]
            if coeff < 0:
                consumes_src.append(mi)
                consumes_dst.append(ri)
                consumes_flux.append(flux_val * abs(coeff))  # flux × stoich
            else:
                produces_src.append(mi)
                produces_dst.append(ri)
                produces_flux.append(flux_val * coeff)

    # GPR edges with gene activity
    gpr_src, gpr_dst, gpr_activity = [], [], []
    for rxn in model.reactions:
        ri = rxn_idx[rxn.id]
        flux_val = abs(fluxes[rxn.id]) if rxn.id in fluxes.index else 0.0
        for gene in rxn.genes:
            gi = gene_idx[gene.id]
            gpr_src.append(gi)
            gpr_dst.append(ri)
            gpr_activity.append(flux_val)

    # Node features: [degree, 1] — same as before
    met_degree = np.zeros(n_met, dtype=np.float32)
    rxn_degree = np.zeros(n_rxn, dtype=np.float32)
    gene_degree = np.zeros(n_gene, dtype=np.float32)
    for rxn in model.reactions:
        ri = rxn_idx[rxn.id]
        rxn_degree[ri] = len(rxn.metabolites) + len(rxn.genes)
        for met in rxn.metabolites:
            met_degree[met_idx[met.id]] += 1
    for gene in model.genes:
        gene_degree[gene_idx[gene.id]] = len(gene.reactions)

    data["metabolite"].x = torch.tensor(
        np.column_stack([met_degree, np.ones(n_met)]), dtype=torch.float32)
    data["reaction"].x = torch.tensor(
        np.column_stack([rxn_degree, np.ones(n_rxn)]), dtype=torch.float32)
    data["gene"].x = torch.tensor(
        np.column_stack([gene_degree, np.ones(n_gene)]), dtype=torch.float32)

    # Edges WITH flux weights (C6=T)
    data["metabolite", "consumes", "reaction"].edge_index = torch.tensor(
        [consumes_src, consumes_dst], dtype=torch.long)
    data["metabolite", "consumes", "reaction"].edge_attr = torch.tensor(
        consumes_flux, dtype=torch.float32).unsqueeze(1)

    data["metabolite", "produces", "reaction"].edge_index = torch.tensor(
        [produces_src, produces_dst], dtype=torch.long)
    data["metabolite", "produces", "reaction"].edge_attr = torch.tensor(
        produces_flux, dtype=torch.float32).unsqueeze(1)

    data["gene", "regulates", "reaction"].edge_index = torch.tensor(
        [gpr_src, gpr_dst], dtype=torch.long)
    data["gene", "regulates", "reaction"].edge_attr = torch.tensor(
        gpr_activity, dtype=torch.float32).unsqueeze(1)

    # Reverse edges
    data["reaction", "rev_consumes", "metabolite"].edge_index = torch.tensor(
        [consumes_dst, consumes_src], dtype=torch.long)
    data["reaction", "rev_consumes", "metabolite"].edge_attr = torch.tensor(
        consumes_flux, dtype=torch.float32).unsqueeze(1)

    data["reaction", "rev_produces", "metabolite"].edge_index = torch.tensor(
        [produces_dst, produces_src], dtype=torch.long)
    data["reaction", "rev_produces", "metabolite"].edge_attr = torch.tensor(
        produces_flux, dtype=torch.float32).unsqueeze(1)

    data["reaction", "rev_regulates", "gene"].edge_index = torch.tensor(
        [gpr_dst, gpr_src], dtype=torch.long)
    data["reaction", "rev_regulates", "gene"].edge_attr = torch.tensor(
        gpr_activity, dtype=torch.float32).unsqueeze(1)

    return data, growth


def build_fixed_graph(model):
    """Build fixed graph (no flux weights) — same for all samples (C2=F, C6=F)."""
    data = HeteroData()

    met_ids = sorted([m.id for m in model.metabolites])
    rxn_ids = sorted([r.id for r in model.reactions])
    gene_ids = sorted([g.id for g in model.genes])

    met_idx = {mid: i for i, mid in enumerate(met_ids)}
    rxn_idx = {rid: i for i, rid in enumerate(rxn_ids)}
    gene_idx = {gid: i for i, gid in enumerate(gene_ids)}

    n_met, n_rxn, n_gene = len(met_ids), len(rxn_ids), len(gene_ids)

    consumes_src, consumes_dst = [], []
    produces_src, produces_dst = [], []
    for rxn in model.reactions:
        ri = rxn_idx[rxn.id]
        for met, coeff in rxn.metabolites.items():
            mi = met_idx[met.id]
            if coeff < 0:
                consumes_src.append(mi); consumes_dst.append(ri)
            else:
                produces_src.append(mi); produces_dst.append(ri)

    gpr_src, gpr_dst = [], []
    for rxn in model.reactions:
        ri = rxn_idx[rxn.id]
        for gene in rxn.genes:
            gpr_src.append(gene_idx[gene.id]); gpr_dst.append(ri)

    met_degree = np.zeros(n_met, dtype=np.float32)
    rxn_degree = np.zeros(n_rxn, dtype=np.float32)
    gene_degree = np.zeros(n_gene, dtype=np.float32)
    for rxn in model.reactions:
        ri = rxn_idx[rxn.id]
        rxn_degree[ri] = len(rxn.metabolites) + len(rxn.genes)
        for met in rxn.metabolites:
            met_degree[met_idx[met.id]] += 1
    for gene in model.genes:
        gene_degree[gene_idx[gene.id]] = len(gene.reactions)

    data["metabolite"].x = torch.tensor(
        np.column_stack([met_degree, np.ones(n_met)]), dtype=torch.float32)
    data["reaction"].x = torch.tensor(
        np.column_stack([rxn_degree, np.ones(n_rxn)]), dtype=torch.float32)
    data["gene"].x = torch.tensor(
        np.column_stack([gene_degree, np.ones(n_gene)]), dtype=torch.float32)

    data["metabolite", "consumes", "reaction"].edge_index = torch.tensor(
        [consumes_src, consumes_dst], dtype=torch.long)
    data["metabolite", "produces", "reaction"].edge_index = torch.tensor(
        [produces_src, produces_dst], dtype=torch.long)
    data["gene", "regulates", "reaction"].edge_index = torch.tensor(
        [gpr_src, gpr_dst], dtype=torch.long)

    data["reaction", "rev_consumes", "metabolite"].edge_index = torch.tensor(
        [consumes_dst, consumes_src], dtype=torch.long)
    data["reaction", "rev_produces", "metabolite"].edge_index = torch.tensor(
        [produces_dst, produces_src], dtype=torch.long)
    data["reaction", "rev_regulates", "gene"].edge_index = torch.tensor(
        [gpr_dst, gpr_src], dtype=torch.long)

    return data


# ============================================================
# 3. HGTGNN with edge attribute support
# ============================================================

class HGTGNNEdgeWeighted(nn.Module):
    """HGTConv that ignores edge_attr (HGTConv doesn't use it natively).
    We encode flux info in node features instead."""

    def __init__(self, metadata, hidden_channels=32, out_channels=32,
                 num_heads=2, num_layers=2, dropout=0.1, node_feature_dim=2):
        super().__init__()
        self.node_types = metadata[0]
        self.edge_types = metadata[1]
        self.out_channels = out_channels
        self.num_layers = num_layers

        self.node_embeddings = nn.ModuleDict()
        for ntype in self.node_types:
            self.node_embeddings[ntype] = nn.Linear(node_feature_dim, hidden_channels)

        self.convs = nn.ModuleList()
        for i in range(num_layers):
            self.convs.append(HGTConv(
                in_channels=hidden_channels,
                out_channels=hidden_channels if i < num_layers - 1 else out_channels,
                metadata=metadata,
                heads=num_heads,
            ))
        self.dropout_layer = nn.Dropout(dropout)

    def forward(self, x_dict, edge_index_dict, knockout_mask=None):
        h_dict = {}
        for ntype in self.node_types:
            if ntype in x_dict:
                feat = x_dict[ntype]
                if ntype == "gene" and knockout_mask is not None:
                    scale = (1.0 - knockout_mask).unsqueeze(1)
                    feat = feat * scale
                h_dict[ntype] = self.node_embeddings[ntype](feat)

        for i, conv in enumerate(self.convs):
            h_prev = {k: v.clone() for k, v in h_dict.items()}
            out_dict = conv(h_dict, edge_index_dict)
            h_dict = {}
            for ntype in self.node_types:
                if ntype in out_dict and out_dict[ntype] is not None:
                    h_dict[ntype] = out_dict[ntype]
                elif ntype in h_prev:
                    h_dict[ntype] = h_prev[ntype]
            if i < self.num_layers - 1:
                h_dict = {k: self.dropout_layer(F.relu(v)) for k, v in h_dict.items()}

        if "gene" in h_dict and h_dict["gene"] is not None:
            graph_emb = h_dict["gene"].mean(dim=0, keepdim=True)
        else:
            embs = [v.mean(dim=0, keepdim=True) for v in h_dict.values() if v is not None]
            graph_emb = torch.cat(embs, dim=1)

        return graph_emb


# ============================================================
# 4. Embedding extraction (MFG version — different graph per sample)
# ============================================================

def extract_mfg_embeddings(gnn, mfg_graphs, knockout_masks):
    """Extract embeddings from condition-specific MFG graphs."""
    gnn.eval()
    all_embeddings = []

    with torch.no_grad():
        for i in range(len(knockout_masks)):
            graph = mfg_graphs[i]
            if graph is None:
                # Infeasible knockout — zero embedding
                all_embeddings.append(np.zeros(gnn.out_channels, dtype=np.float32))
                continue

            x_dict = {
                "metabolite": graph["metabolite"].x,
                "reaction": graph["reaction"].x,
                "gene": graph["gene"].x,
            }
            edge_index_dict = {etype: graph[etype].edge_index for etype in graph.edge_types}
            ko_mask = torch.tensor(knockout_masks[i], dtype=torch.float32)

            emb = gnn(x_dict, edge_index_dict, knockout_mask=ko_mask)
            all_embeddings.append(emb.squeeze(0).numpy())

    return np.array(all_embeddings)


def extract_fixed_embeddings(gnn, graph, knockout_masks):
    """Extract embeddings from fixed graph (same for all samples)."""
    gnn.eval()
    x_dict = {
        "metabolite": graph["metabolite"].x,
        "reaction": graph["reaction"].x,
        "gene": graph["gene"].x,
    }
    edge_index_dict = {etype: graph[etype].edge_index for etype in graph.edge_types}

    all_embeddings = []
    with torch.no_grad():
        for i in range(len(knockout_masks)):
            ko_mask = torch.tensor(knockout_masks[i], dtype=torch.float32)
            emb = gnn(x_dict, edge_index_dict, knockout_mask=ko_mask)
            all_embeddings.append(emb.squeeze(0).numpy())

    return np.array(all_embeddings)


# ============================================================
# 5. Main Experiment
# ============================================================

def run_flowgat_experiment():
    print("=" * 70)
    print("NAP FlowGAT-style Experiment: C2+C6 Validation")
    print("=" * 70)
    print()
    print("NAP Predictions:")
    print("  A. XGBoost-only:        C2=F, C6=F → 0/6 → GNN no value expected")
    print("  B. GNN+XGB (fixed):     C2=F, C6=F → 0/6 → GNN no value expected")
    print("  C. GNN+XGB (MFG):       C2=T, C6=T → 2/6 → GNN value EXPECTED")
    print()

    model = load_textbook()
    wt_growth = model.optimize().objective_value
    print(f"WT growth: {wt_growth:.4f}")

    # --- Generate data ---
    N_SAMPLES = 150  # Keep small for speed
    print(f"\n[1] Generating {N_SAMPLES} knockout combos + FBA + MFG...")

    combos = generate_random_combos(model, N_SAMPLES, min_k=1, max_k=5, seed=42)
    knockout_masks = build_knockout_mask(model, combos)

    # Run FBA and build MFG for each knockout
    growth_rates = []
    mfg_graphs = []
    t0 = time.time()

    for i, combo in enumerate(combos):
        mfg, growth = build_mfg_for_knockout(model, [g.id for g in model.genes], combo)
        growth_rates.append(growth)
        mfg_graphs.append(mfg)

        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{N_SAMPLES} done ({time.time()-t0:.1f}s)")

    growth_rates = np.array(growth_rates, dtype=np.float32)
    n_feasible = (growth_rates > 1e-6).sum()
    print(f"  Done: {n_feasible}/{N_SAMPLES} feasible, {time.time()-t0:.1f}s")

    # Build fixed graph
    fixed_graph = build_fixed_graph(model)

    # --- Experiment A: XGBoost-only ---
    print("\n[2] Experiment A: XGBoost-only baseline (5-fold CV)...")

    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    xgb_r2_scores = []

    for fold, (train_idx, test_idx) in enumerate(kf.split(knockout_masks)):
        X_train, X_test = knockout_masks[train_idx], knockout_masks[test_idx]
        y_train, y_test = growth_rates[train_idx], growth_rates[test_idx]

        xgb_model = xgb.XGBRegressor(
            objective='reg:squarederror', max_depth=6,
            learning_rate=0.1, n_estimators=200,
            subsample=0.8, colsample_bytree=0.8,
            random_state=42, verbosity=0,
        )
        xgb_model.fit(X_train, y_train)
        y_pred = xgb_model.predict(X_test)
        xgb_r2_scores.append(r2_score(y_test, y_pred))

    xgb_r2 = np.mean(xgb_r2_scores)
    xgb_r2_std = np.std(xgb_r2_scores)
    print(f"  XGBoost-only: R2 = {xgb_r2:.4f} +/- {xgb_r2_std:.4f}")

    # --- Experiment B: GNN+XGBoost (fixed graph) ---
    print("\n[3] Experiment B: GNN+XGBoost (fixed graph, C2=F, C6=F)...")

    metadata = fixed_graph.metadata()
    gnn_fixed = HGTGNNEdgeWeighted(metadata, hidden_channels=32, out_channels=32,
                                    num_heads=2, num_layers=2)

    # Train GNN supervised on growth rate
    x_dict = {
        "metabolite": fixed_graph["metabolite"].x,
        "reaction": fixed_graph["reaction"].x,
        "gene": fixed_graph["gene"].x,
    }
    edge_index_dict = {etype: fixed_graph[etype].edge_index for etype in fixed_graph.edge_types}
    reg_head = nn.Linear(32, 1)
    optimizer = torch.optim.Adam(list(gnn_fixed.parameters()) + list(reg_head.parameters()), lr=0.01)

    masks_t = torch.tensor(knockout_masks, dtype=torch.float32)
    y_t = torch.tensor(growth_rates, dtype=torch.float32)

    gnn_fixed.train()
    for epoch in range(30):
        optimizer.zero_grad()
        total_loss = 0
        for i in range(len(masks_t)):
            emb = gnn_fixed(x_dict, edge_index_dict, knockout_mask=masks_t[i])
            pred = reg_head(emb).squeeze()
            loss = F.mse_loss(pred, y_t[i])
            loss.backward(retain_graph=True)
            total_loss += loss.item()
        optimizer.step()

    # Extract fixed-graph embeddings and evaluate
    fixed_embs = extract_fixed_embeddings(gnn_fixed, fixed_graph, knockout_masks)
    combined_fixed = np.concatenate([fixed_embs, knockout_masks], axis=1)

    fixed_r2_scores = []
    for fold, (train_idx, test_idx) in enumerate(kf.split(combined_fixed)):
        X_train, X_test = combined_fixed[train_idx], combined_fixed[test_idx]
        y_train, y_test = growth_rates[train_idx], growth_rates[test_idx]

        xgb_model = xgb.XGBRegressor(
            objective='reg:squarederror', max_depth=6,
            learning_rate=0.1, n_estimators=200,
            random_state=42, verbosity=0,
        )
        xgb_model.fit(X_train, y_train)
        y_pred = xgb_model.predict(X_test)
        fixed_r2_scores.append(r2_score(y_test, y_pred))

    fixed_r2 = np.mean(fixed_r2_scores)
    fixed_r2_std = np.std(fixed_r2_scores)
    print(f"  GNN(fixed)+XGB: R2 = {fixed_r2:.4f} +/- {fixed_r2_std:.4f}")

    # --- Experiment C: GNN+XGBoost (MFG, C2=T, C6=T) ---
    print("\n[4] Experiment C: GNN+XGBoost (MFG, C2=T, C6=T)...")
    print("  Note: Each sample has its own graph (flux-weighted edges)")
    print("  Training GNN on MFG-conditioned samples...")

    # For MFG: train GNN on a few MFG graphs, extract embeddings for all
    gnn_mfg = HGTGNNEdgeWeighted(metadata, hidden_channels=32, out_channels=32,
                                  num_heads=2, num_layers=2)
    reg_head_mfg = nn.Linear(32, 1)
    optimizer_mfg = torch.optim.Adam(list(gnn_mfg.parameters()) + list(reg_head_mfg.parameters()), lr=0.01)

    gnn_mfg.train()
    for epoch in range(30):
        optimizer_mfg.zero_grad()
        total_loss = 0
        n_valid = 0
        for i in range(len(knockout_masks)):
            if mfg_graphs[i] is None:
                continue

            g = mfg_graphs[i]
            x_dict_mfg = {
                "metabolite": g["metabolite"].x,
                "reaction": g["reaction"].x,
                "gene": g["gene"].x,
            }
            edge_dict_mfg = {etype: g[etype].edge_index for etype in g.edge_types}
            ko_mask = masks_t[i]

            emb = gnn_mfg(x_dict_mfg, edge_dict_mfg, knockout_mask=ko_mask)
            pred = reg_head_mfg(emb).squeeze()
            loss = F.mse_loss(pred, y_t[i])
            loss.backward(retain_graph=True)
            total_loss += loss.item()
            n_valid += 1
        optimizer_mfg.step()

        if (epoch + 1) % 10 == 0:
            print(f"    Epoch {epoch+1}: avg_loss={total_loss/max(n_valid,1):.4f} ({n_valid} valid)")

    # Extract MFG embeddings
    print("  Extracting MFG embeddings...")
    mfg_embs = extract_mfg_embeddings(gnn_mfg, mfg_graphs, knockout_masks)
    combined_mfg = np.concatenate([mfg_embs, knockout_masks], axis=1)

    mfg_r2_scores = []
    for fold, (train_idx, test_idx) in enumerate(kf.split(combined_mfg)):
        X_train, X_test = combined_mfg[train_idx], combined_mfg[test_idx]
        y_train, y_test = growth_rates[train_idx], growth_rates[test_idx]

        xgb_model = xgb.XGBRegressor(
            objective='reg:squarederror', max_depth=6,
            learning_rate=0.1, n_estimators=200,
            random_state=42, verbosity=0,
        )
        xgb_model.fit(X_train, y_train)
        y_pred = xgb_model.predict(X_test)
        mfg_r2_scores.append(r2_score(y_test, y_pred))

    mfg_r2 = np.mean(mfg_r2_scores)
    mfg_r2_std = np.std(mfg_r2_scores)
    print(f"  GNN(MFG)+XGB: R2 = {mfg_r2:.4f} +/- {mfg_r2_std:.4f}")

    # --- Experiment D: Flux features as tabular (ablation) ---
    # Test: what if we just add flux summary statistics as tabular features?
    print("\n[5] Experiment D: Flux statistics as tabular features (ablation)...")

    # For each knockout, compute flux summary stats
    flux_features = []
    for i, mfg in enumerate(mfg_graphs):
        if mfg is None:
            flux_features.append(np.zeros(5, dtype=np.float32))
            continue

        # Extract flux-weighted edge statistics
        consumes_flux = mfg["metabolite", "consumes", "reaction"].edge_attr.squeeze().numpy()
        produces_flux = mfg["metabolite", "produces", "reaction"].edge_attr.squeeze().numpy()
        gpr_flux = mfg["gene", "regulates", "reaction"].edge_attr.squeeze().numpy()

        all_flux = np.concatenate([consumes_flux, produces_flux, gpr_flux])
        flux_features.append([
            all_flux.mean(),
            all_flux.std(),
            all_flux.max(),
            all_flux.min(),
            (all_flux > 0.01).sum() / len(all_flux),  # active edge fraction
        ])

    flux_features = np.array(flux_features, dtype=np.float32)
    combined_flux = np.concatenate([knockout_masks, flux_features], axis=1)

    flux_r2_scores = []
    for fold, (train_idx, test_idx) in enumerate(kf.split(combined_flux)):
        X_train, X_test = combined_flux[train_idx], combined_flux[test_idx]
        y_train, y_test = growth_rates[train_idx], growth_rates[test_idx]

        xgb_model = xgb.XGBRegressor(
            objective='reg:squarederror', max_depth=6,
            learning_rate=0.1, n_estimators=200,
            random_state=42, verbosity=0,
        )
        xgb_model.fit(X_train, y_train)
        y_pred = xgb_model.predict(X_test)
        flux_r2_scores.append(r2_score(y_test, y_pred))

    flux_r2 = np.mean(flux_r2_scores)
    flux_r2_std = np.std(flux_r2_scores)
    print(f"  XGB+FluxStats: R2 = {flux_r2:.4f} +/- {flux_r2_std:.4f}")

    # --- Summary ---
    print("\n" + "=" * 70)
    print("NAP PREDICTION vs EXPERIMENTAL RESULT")
    print("=" * 70)
    print(f"{'Experiment':<30} {'NAP':<12} {'R2':<12} {'R2_std':<12} {'vs XGB':<12}")
    print("-" * 78)
    print(f"{'A. XGBoost-only':<30} {'0/6':<12} {xgb_r2:<12.4f} {xgb_r2_std:<12.4f} {'baseline':<12}")
    print(f"{'B. GNN(fixed)+XGB':<30} {'0/6':<12} {fixed_r2:<12.4f} {fixed_r2_std:<12.4f} {fixed_r2-xgb_r2:<+12.4f}")
    print(f"{'C. GNN(MFG)+XGB':<30} {'2/6':<12} {mfg_r2:<12.4f} {mfg_r2_std:<12.4f} {mfg_r2-xgb_r2:<+12.4f}")
    print(f"{'D. XGB+FluxStats':<30} {'ablation':<12} {flux_r2:<12.4f} {flux_r2_std:<12.4f} {flux_r2-xgb_r2:<+12.4f}")
    print()

    # NAP verdict
    print("NAP VERDICT:")
    if mfg_r2 > xgb_r2 + 0.02:
        print("  C2+C6 conditions CONFIRMED: MFG-based GNN adds value over XGBoost-only")
        print("  NAP prediction (2/6 → GNN value) is CORRECT")
    elif mfg_r2 > xgb_r2:
        print("  C2+C6 conditions MARGINAL: MFG-based GNN slightly better")
        print("  NAP prediction directionally correct but effect size small")
    else:
        print("  C2+C6 conditions NOT CONFIRMED: MFG-based GNN does NOT add value")
        print("  NAP prediction may be INCORRECT or experiment underpowered")

    if flux_r2 > xgb_r2 + 0.02:
        print("  Flux statistics as tabular features also add value")
        print("  → GNN may not be necessary; engineered features may suffice")

    print()
    print("CIRCULAR LOGIC NOTE:")
    print("  MFG construction requires FBA, but surrogate's purpose is to replace FBA.")
    print("  Practical solution: coarse FBA/pFBA for MFG + surrogate for growth.")
    print("  This experiment tests theoretical value, not practical deployability.")

    return {
        "xgb_r2": xgb_r2,
        "fixed_r2": fixed_r2,
        "mfg_r2": mfg_r2,
        "flux_r2": flux_r2,
    }


if __name__ == "__main__":
    results = run_flowgat_experiment()
