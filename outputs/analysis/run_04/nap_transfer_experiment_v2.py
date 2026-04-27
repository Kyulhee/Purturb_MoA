"""
NAP Transfer Experiment v2: Common Embedding Space for Cross-Species GNN Transfer
==================================================================================
Key Fix (from v1): Cross-species transfer requires mapping to a COMMON embedding
space, not direct HGTConv weight transfer.

Architecture:
  Source species: SourceGNN(graph_s) → 32d → SharedProjection → common_dim → RegHead → growth
  Target species: TargetGNN(graph_t) → 32d → SharedProjection → common_dim → RegHead → growth

  Transfer = freeze SharedProjection + RegHead from source, only train TargetGNN encoder.
  This tests whether the "common metabolic embedding space" learned from source
  provides a better starting point for the target species.

NAP Prediction:
  - textbook alone: C1=F, C2=F → 0/6 → GNN no value
  - textbook→iSB619: C2=T(graph differs), C3=T(transfer to new domain) → 2/6 → GNN value
  - textbook→iCN718: C2=T, C3=T → 2/6 → GNN value

Key Question: Does transfer via shared projection improve over XGBoost-only on target?
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
from torch_geometric.nn import HGTConv
import xgboost as xgb
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from scipy.spatial.distance import cdist

warnings.filterwarnings("ignore")

# ============================================================
# 1. Model Loading & FBA Data Generation
# ============================================================

def load_model(model_id: str) -> cobra.Model:
    """Load a BiGG model by ID."""
    model = cobra.io.load_model(model_id)
    print(f"[Model] {model_id}: {len(model.reactions)} rxn, "
          f"{len(model.metabolites)} met, {len(model.genes)} genes")
    return model


def run_fba_knockout(model: cobra.Model, gene_ids: List[str]) -> float:
    """Single FBA knockout, return growth rate."""
    try:
        model_cp = model.copy()
        with model_cp:
            for gid in gene_ids:
                if gid in [g.id for g in model_cp.genes]:
                    model_cp.genes.get_by_id(gid).knock_out()
            sol = model_cp.optimize()
            if sol.status == "optimal":
                return sol.objective_value
        return 0.0
    except:
        return 0.0


def generate_random_combos(model: cobra.Model, n: int,
                           min_k: int = 1, max_k: int = 5,
                           seed: int = 42) -> List[List[str]]:
    """Generate random knockout combinations."""
    gene_ids = [g.id for g in model.genes]
    rng = random.Random(seed)
    combos = []
    for _ in range(n):
        k = rng.randint(min_k, min(max_k, len(gene_ids)))
        combo = rng.sample(gene_ids, k)
        combos.append(combo)
    return combos


def build_knockout_mask(model: cobra.Model, combos: List[List[str]]) -> np.ndarray:
    """Build binary knockout mask (n_combos, n_genes)."""
    gene_ids = sorted([g.id for g in model.genes])
    gene_idx = {gid: i for i, gid in enumerate(gene_ids)}
    n_gene = len(gene_ids)
    mask = np.zeros((len(combos), n_gene), dtype=np.float32)
    for i, combo in enumerate(combos):
        for gid in combo:
            if gid in gene_idx:
                mask[i, gene_idx[gid]] = 1.0
    return mask


def generate_fba_dataset(model: cobra.Model, n_samples: int = 500,
                         seed: int = 42) -> Tuple[np.ndarray, np.ndarray]:
    """Generate FBA knockout dataset for a model."""
    t0 = time.time()
    combos = generate_random_combos(model, n_samples, seed=seed)
    masks = build_knockout_mask(model, combos)
    growths = np.array([run_fba_knockout(model, c) for c in combos],
                       dtype=np.float32)
    elapsed = time.time() - t0
    n_feasible = (growths > 1e-6).sum()
    print(f"  FBA dataset: {n_samples} samples, {n_feasible} feasible, "
          f"{elapsed:.1f}s ({elapsed/n_samples*1000:.0f}ms/call)")
    return masks, growths


# ============================================================
# 2. Graph Conversion
# ============================================================

def model_to_heterodata(model: cobra.Model) -> HeteroData:
    """Convert COBRApy model to PyG HeteroData."""
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
                consumes_src.append(mi)
                consumes_dst.append(ri)
            else:
                produces_src.append(mi)
                produces_dst.append(ri)

    gpr_src, gpr_dst = [], []
    for rxn in model.reactions:
        ri = rxn_idx[rxn.id]
        for gene in rxn.genes:
            gi = gene_idx[gene.id]
            gpr_src.append(gi)
            gpr_dst.append(ri)

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
# 3. GNN Encoder + Shared Projection Head
# ============================================================

class HGTGNNEncoder(nn.Module):
    """Species-specific HGT encoder: graph → per-sample 32d embedding."""

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

        # Graph-level readout: mean-pool gene embeddings
        if "gene" in h_dict and h_dict["gene"] is not None:
            graph_emb = h_dict["gene"].mean(dim=0, keepdim=True)
        else:
            embs = [v.mean(dim=0, keepdim=True) for v in h_dict.values() if v is not None]
            graph_emb = torch.cat(embs, dim=1)

        return graph_emb  # (1, out_channels)


class SharedProjectionHead(nn.Module):
    """
    Shared projection that maps species-specific GNN embeddings to a COMMON
    embedding space. This is the key component that enables cross-species transfer.

    Architecture: gnn_dim → common_dim → common_dim (with ReLU + LayerNorm)
    """

    def __init__(self, gnn_dim: int = 32, common_dim: int = 16):
        super().__init__()
        self.proj = nn.Sequential(
            nn.Linear(gnn_dim, common_dim),
            nn.ReLU(),
            nn.LayerNorm(common_dim),
            nn.Linear(common_dim, common_dim),
        )
        self.common_dim = common_dim

    def forward(self, gnn_embedding):
        """
        Args:
            gnn_embedding: (batch_size, gnn_dim) or (1, gnn_dim)
        Returns:
            common_embedding: (batch_size, common_dim)
        """
        return self.proj(gnn_embedding)


class GrowthRegressionHead(nn.Module):
    """Regression head from common embedding to growth rate."""

    def __init__(self, common_dim: int = 16):
        super().__init__()
        self.reg = nn.Sequential(
            nn.Linear(common_dim, 8),
            nn.ReLU(),
            nn.Linear(8, 1),
        )

    def forward(self, common_embedding):
        return self.reg(common_embedding).squeeze(-1)


class MetabolicSurrogateModel(nn.Module):
    """
    Full model: GNN encoder (species-specific) + SharedProjection + RegHead.

    For SOURCE training: all three components are trained end-to-end.
    For TRANSFER: SharedProjection + RegHead are frozen from source,
                  only the target GNN encoder is trained.
    """

    def __init__(self, encoder: HGTGNNEncoder,
                 projection: SharedProjectionHead,
                 regression: GrowthRegressionHead):
        super().__init__()
        self.encoder = encoder
        self.projection = projection
        self.regression = regression

    def forward(self, x_dict, edge_index_dict, knockout_mask=None):
        gnn_emb = self.encoder(x_dict, edge_index_dict, knockout_mask)
        common_emb = self.projection(gnn_emb)
        growth_pred = self.regression(common_emb)
        return growth_pred, common_emb


# ============================================================
# 4. Training Utilities
# ============================================================

def train_end_to_end(model: MetabolicSurrogateModel,
                     graph: HeteroData,
                     masks: np.ndarray,
                     growth: np.ndarray,
                     epochs: int = 30,
                     lr: float = 0.005,
                     batch_size: int = 32,
                     verbose: bool = True) -> MetabolicSurrogateModel:
    """Train full model end-to-end on a single species."""
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    x_dict = {
        "metabolite": graph["metabolite"].x,
        "reaction": graph["reaction"].x,
        "gene": graph["gene"].x,
    }
    edge_index_dict = {etype: graph[etype].edge_index for etype in graph.edge_types}

    masks_t = torch.tensor(masks, dtype=torch.float32)
    y_t = torch.tensor(growth, dtype=torch.float32)
    n_samples = len(masks_t)

    model.train()
    for epoch in range(epochs):
        total_loss = 0.0
        n_batches = 0
        perm = torch.randperm(n_samples)

        for start in range(0, n_samples, batch_size):
            idx = perm[start:start + batch_size]
            optimizer.zero_grad()

            batch_loss = 0.0
            for i in idx:
                pred, _ = model(x_dict, edge_index_dict, knockout_mask=masks_t[i])
                batch_loss += F.mse_loss(pred, y_t[i])

            batch_loss = batch_loss / len(idx)
            batch_loss.backward()
            # Gradient clipping
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            total_loss += batch_loss.item()
            n_batches += 1

        scheduler.step()

        if verbose and (epoch + 1) % 10 == 0:
            print(f"    Epoch {epoch+1}/{epochs}: loss={total_loss/n_batches:.4f}")

    return model


def train_encoder_only(model: MetabolicSurrogateModel,
                       graph: HeteroData,
                       masks: np.ndarray,
                       growth: np.ndarray,
                       epochs: int = 30,
                       lr: float = 0.005,
                       batch_size: int = 32,
                       verbose: bool = True) -> MetabolicSurrogateModel:
    """
    Train ONLY the GNN encoder, with SharedProjection + RegHead FROZEN.
    This simulates transfer: source-trained projection/regression is reused,
    only the target-specific encoder adapts.
    """
    # Freeze projection and regression heads
    for param in model.projection.parameters():
        param.requires_grad = False
    for param in model.regression.parameters():
        param.requires_grad = False

    optimizer = torch.optim.Adam(model.encoder.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    x_dict = {
        "metabolite": graph["metabolite"].x,
        "reaction": graph["reaction"].x,
        "gene": graph["gene"].x,
    }
    edge_index_dict = {etype: graph[etype].edge_index for etype in graph.edge_types}

    masks_t = torch.tensor(masks, dtype=torch.float32)
    y_t = torch.tensor(growth, dtype=torch.float32)
    n_samples = len(masks_t)

    model.train()
    # Even though projection/regression are frozen, we need them in eval mode
    # for stable outputs during encoder training
    model.projection.eval()
    model.regression.eval()

    for epoch in range(epochs):
        total_loss = 0.0
        n_batches = 0
        perm = torch.randperm(n_samples)

        for start in range(0, n_samples, batch_size):
            idx = perm[start:start + batch_size]
            optimizer.zero_grad()

            batch_loss = 0.0
            for i in idx:
                pred, _ = model(x_dict, edge_index_dict, knockout_mask=masks_t[i])
                batch_loss += F.mse_loss(pred, y_t[i])

            batch_loss = batch_loss / len(idx)
            batch_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.encoder.parameters(), max_norm=1.0)
            optimizer.step()
            total_loss += batch_loss.item()
            n_batches += 1

        scheduler.step()

        if verbose and (epoch + 1) % 10 == 0:
            print(f"    Epoch {epoch+1}/{epochs}: loss={total_loss/n_batches:.4f}")

    # Unfreeze for later use
    for param in model.projection.parameters():
        param.requires_grad = True
    for param in model.regression.parameters():
        param.requires_grad = True

    return model


def extract_common_embeddings(model: MetabolicSurrogateModel,
                              graph: HeteroData,
                              masks: np.ndarray,
                              batch_size: int = 64) -> np.ndarray:
    """Extract common-space embeddings for a batch of knockout masks."""
    model.eval()
    x_dict = {
        "metabolite": graph["metabolite"].x,
        "reaction": graph["reaction"].x,
        "gene": graph["gene"].x,
    }
    edge_index_dict = {etype: graph[etype].edge_index for etype in graph.edge_types}

    masks_t = torch.tensor(masks, dtype=torch.float32)
    all_embeddings = []

    with torch.no_grad():
        for i in range(0, len(masks_t), batch_size):
            batch_masks = masks_t[i:i + batch_size]
            for j in range(len(batch_masks)):
                _, common_emb = model(x_dict, edge_index_dict,
                                      knockout_mask=batch_masks[j])
                all_embeddings.append(common_emb.squeeze(0).numpy())

    return np.array(all_embeddings)


# ============================================================
# 5. Experiment Functions
# ============================================================

def experiment_xgboost_only(masks, growth, n_folds=5):
    """XGBoost-only baseline with k-fold CV."""
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
    r2_scores = []

    for train_idx, test_idx in kf.split(masks):
        X_train, X_test = masks[train_idx], masks[test_idx]
        y_train, y_test = growth[train_idx], growth[test_idx]

        model = xgb.XGBRegressor(
            objective='reg:squarederror', max_depth=6,
            learning_rate=0.1, n_estimators=200,
            subsample=0.8, colsample_bytree=0.8,
            random_state=42, verbosity=0,
        )
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        r2_scores.append(r2_score(y_test, y_pred))

    return np.mean(r2_scores), np.std(r2_scores)


def experiment_gnn_end_to_end(graph, masks, growth, n_folds=3,
                               hidden=32, out=32, common_dim=16,
                               epochs=30, lr=0.005):
    """GNN+SharedProjection+RegHead trained end-to-end on a single species."""
    metadata = graph.metadata()
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
    r2_scores = []

    for fold, (train_idx, test_idx) in enumerate(kf.split(masks)):
        X_train, X_test = masks[train_idx], masks[test_idx]
        y_train, y_test = growth[train_idx], growth[test_idx]

        encoder = HGTGNNEncoder(metadata, hidden_channels=hidden,
                                 out_channels=out, num_heads=2, num_layers=2)
        projection = SharedProjectionHead(gnn_dim=out, common_dim=common_dim)
        regression = GrowthRegressionHead(common_dim=common_dim)
        model = MetabolicSurrogateModel(encoder, projection, regression)

        model = train_end_to_end(model, graph, X_train, y_train,
                                  epochs=epochs, lr=lr, verbose=True)

        # Evaluate
        model.eval()
        x_dict = {
            "metabolite": graph["metabolite"].x,
            "reaction": graph["reaction"].x,
            "gene": graph["gene"].x,
        }
        edge_index_dict = {etype: graph[etype].edge_index for etype in graph.edge_types}
        masks_test = torch.tensor(X_test, dtype=torch.float32)

        preds = []
        with torch.no_grad():
            for i in range(len(masks_test)):
                pred, _ = model(x_dict, edge_index_dict, knockout_mask=masks_test[i])
                preds.append(pred.item())

        r2 = r2_score(y_test, np.array(preds))
        r2_scores.append(r2)
        print(f"    Fold {fold+1}: End-to-End R2 = {r2:.4f}")

    return np.mean(r2_scores), np.std(r2_scores)


def experiment_transfer(source_model: MetabolicSurrogateModel,
                        target_graph: HeteroData,
                        masks_target: np.ndarray,
                        growth_target: np.ndarray,
                        n_folds=3,
                        fine_tune_epochs=30,
                        lr=0.005,
                        common_dim=16):
    """
    TRANSFER via common embedding space:
    1. Source model's SharedProjection + RegHead are TRANSFERRED (frozen initially)
    2. Target-specific GNN encoder is trained to produce embeddings
       that are meaningful in the source-learned common space.
    3. Then optionally unfreeze projection for fine-tuning.

    Compare against:
    - Scratch: same architecture but all trained from scratch on target
    - XGBoost-only: tabular baseline
    """
    metadata = target_graph.metadata()

    kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
    r2_transfer = []
    r2_scratch = []

    for fold, (train_idx, test_idx) in enumerate(kf.split(masks_target)):
        X_train, X_test = masks_target[train_idx], masks_target[test_idx]
        y_train, y_test = growth_target[train_idx], growth_target[test_idx]

        # === (A) FROM SCRATCH: all components trained from scratch on target ===
        encoder_scratch = HGTGNNEncoder(metadata, hidden_channels=32,
                                         out_channels=32, num_heads=2, num_layers=2)
        projection_scratch = SharedProjectionHead(gnn_dim=32, common_dim=common_dim)
        regression_scratch = GrowthRegressionHead(common_dim=common_dim)
        model_scratch = MetabolicSurrogateModel(
            encoder_scratch, projection_scratch, regression_scratch)

        model_scratch = train_end_to_end(
            model_scratch, target_graph, X_train, y_train,
            epochs=fine_tune_epochs, lr=lr, verbose=False)

        # Evaluate scratch
        model_scratch.eval()
        x_dict = {
            "metabolite": target_graph["metabolite"].x,
            "reaction": target_graph["reaction"].x,
            "gene": target_graph["gene"].x,
        }
        edge_index_dict = {etype: target_graph[etype].edge_index
                           for etype in target_graph.edge_types}
        masks_test = torch.tensor(X_test, dtype=torch.float32)

        preds_s = []
        with torch.no_grad():
            for i in range(len(masks_test)):
                pred, _ = model_scratch(x_dict, edge_index_dict,
                                        knockout_mask=masks_test[i])
                preds_s.append(pred.item())
        r2_s = r2_score(y_test, np.array(preds_s))
        r2_scratch.append(r2_s)

        # === (B) TRANSFER: use source projection + regression, train target encoder ===
        # Create new encoder for target species
        encoder_transfer = HGTGNNEncoder(metadata, hidden_channels=32,
                                          out_channels=32, num_heads=2, num_layers=2)

        # Clone source projection + regression (these define the common embedding space)
        projection_transfer = SharedProjectionHead(gnn_dim=32, common_dim=common_dim)
        regression_transfer = GrowthRegressionHead(common_dim=common_dim)

        # COPY source weights to projection + regression
        projection_transfer.load_state_dict(source_model.projection.state_dict())
        regression_transfer.load_state_dict(source_model.regression.state_dict())

        model_transfer = MetabolicSurrogateModel(
            encoder_transfer, projection_transfer, regression_transfer)

        # Phase 1: Train ONLY the encoder (projection/regression frozen from source)
        model_transfer = train_encoder_only(
            model_transfer, target_graph, X_train, y_train,
            epochs=fine_tune_epochs // 2, lr=lr, verbose=False)

        # Phase 2: Fine-tune ALL components with lower learning rate
        model_transfer = train_end_to_end(
            model_transfer, target_graph, X_train, y_train,
            epochs=fine_tune_epochs // 2, lr=lr * 0.1, verbose=False)

        # Evaluate transfer
        model_transfer.eval()
        preds_t = []
        with torch.no_grad():
            for i in range(len(masks_test)):
                pred, _ = model_transfer(x_dict, edge_index_dict,
                                         knockout_mask=masks_test[i])
                preds_t.append(pred.item())
        r2_t = r2_score(y_test, np.array(preds_t))
        r2_transfer.append(r2_t)

        print(f"    Fold {fold+1}: Scratch R2={r2_s:.4f}, Transfer R2={r2_t:.4f}, "
              f"Delta={r2_t - r2_s:+.4f}")

    return (np.mean(r2_scratch), np.std(r2_scratch),
            np.mean(r2_transfer), np.std(r2_transfer))


def experiment_gnn_xgboost_hybrid(graph, masks, growth, n_folds=3,
                                   hidden=32, out=32, common_dim=16,
                                   epochs=30, lr=0.005):
    """GNN common embeddings + XGBoost hybrid (alternative to end-to-end regression)."""
    metadata = graph.metadata()
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
    r2_scores = []

    for fold, (train_idx, test_idx) in enumerate(kf.split(masks)):
        X_train, X_test = masks[train_idx], masks[test_idx]
        y_train, y_test = growth[train_idx], growth[test_idx]

        encoder = HGTGNNEncoder(metadata, hidden_channels=hidden,
                                 out_channels=out, num_heads=2, num_layers=2)
        projection = SharedProjectionHead(gnn_dim=out, common_dim=common_dim)
        regression = GrowthRegressionHead(common_dim=common_dim)
        model = MetabolicSurrogateModel(encoder, projection, regression)

        model = train_end_to_end(model, graph, X_train, y_train,
                                  epochs=epochs, lr=lr, verbose=False)

        # Extract common embeddings
        train_embs = extract_common_embeddings(model, graph, X_train)
        test_embs = extract_common_embeddings(model, graph, X_test)

        # Hybrid: common embeddings + knockout mask → XGBoost
        combined_train = np.concatenate([train_embs, X_train], axis=1)
        combined_test = np.concatenate([test_embs, X_test], axis=1)

        xgb_model = xgb.XGBRegressor(
            objective='reg:squarederror', max_depth=6,
            learning_rate=0.1, n_estimators=200,
            random_state=42, verbosity=0,
        )
        xgb_model.fit(combined_train, y_train)
        y_pred = xgb_model.predict(combined_test)
        r2 = r2_score(y_test, y_pred)
        r2_scores.append(r2)

        # Also test embeddings-only XGBoost (no knockout mask)
        xgb_emb = xgb.XGBRegressor(
            objective='reg:squarederror', max_depth=4,
            learning_rate=0.1, n_estimators=100,
            random_state=42, verbosity=0,
        )
        xgb_emb.fit(train_embs, y_train)
        r2_emb = r2_score(y_test, xgb_emb.predict(test_embs))

        print(f"    Fold {fold+1}: Hybrid R2={r2:.4f}, Emb-only R2={r2_emb:.4f}")

    return np.mean(r2_scores), np.std(r2_scores)


# ============================================================
# 6. Feature Completeness Analysis (C1 condition)
# ============================================================

def analyze_feature_completeness(masks, growth):
    """Check if knockout mask fully determines growth rate (C1 condition)."""
    from sklearn.feature_selection import mutual_info_regression

    mi = mutual_info_regression(masks, growth, random_state=42)
    total_mi = mi.sum()
    h_y = np.var(growth)
    ratio = total_mi / h_y if h_y > 0 else 0
    print(f"  Feature completeness: MI={total_mi:.4f}, Var(Y)={h_y:.4f}, ratio={ratio:.4f}")
    print(f"  C1 (feature incomplete): {'T' if ratio < 0.9 else 'F'} (ratio < 0.9)")
    return ratio


# ============================================================
# 7. Common Embedding Space Alignment Analysis
# ============================================================

def analyze_embedding_alignment(source_model, target_model,
                                 source_graph, target_graph,
                                 source_masks, target_masks,
                                 n_samples=50):
    """
    Analyze whether source and target common embeddings align.
    If they do, transfer should work (C3 validated).
    If they don't, transfer won't help (C3 not sufficient).
    """
    source_model.eval()
    target_model.eval()

    # Extract common embeddings for both species
    source_embs = extract_common_embeddings(
        source_model, source_graph, source_masks[:n_samples])
    target_embs = extract_common_embeddings(
        target_model, target_graph, target_masks[:n_samples])

    # CCA (Canonical Correlation Analysis) to measure alignment
    from sklearn.cross_decomposition import CCA
    n_components = min(source_embs.shape[1], target_embs.shape[1], 8)
    cca = CCA(n_components=n_components)
    source_c, target_c = cca.fit_transform(source_embs, target_embs)

    correlations = []
    for i in range(n_components):
        corr = np.corrcoef(source_c[:, i], target_c[:, i])[0, 1]
        correlations.append(corr)

    mean_corr = np.mean(correlations)
    print(f"  CCA alignment: mean corr = {mean_corr:.4f}, "
          f"per-component = {[f'{c:.3f}' for c in correlations]}")
    print(f"  C3 (transfer viable): {'T' if mean_corr > 0.3 else 'F'} "
          f"(mean CCA corr > 0.3)")

    return mean_corr, correlations
