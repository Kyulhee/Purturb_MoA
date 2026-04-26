"""
Module B: GNN + XGBoost Surrogate Model
=========================================
HGTConv-based heterogeneous GNN encoder + XGBoost regressor.
GNN produces graph-level embeddings from the metabolic network,
which are concatenated with the knockout mask and fed to XGBoost.

Architecture:
  HeteroData → HGTConv(2-layer) → gene-mean-pool → graph_emb(32d)
  graph_emb(32d) ⊕ knockout_mask(137d) → XGBoost → growth_rate

Key design:
  - Knockout mask is NOT added to gene features (avoids dimension mismatch)
  - Instead, knockout mask modifies gene features via element-wise multiply
  - GNN processes the modified graph → embedding
  - Embedding concatenated with raw knockout mask → XGBoost

Design decisions (from stages/03_planning.md):
  - HGTConv: supports 3 node types (metabolite/reaction/gene) natively
  - Bidirectional edges: reverse edges (reaction→metabolite, reaction→gene) are
    required so that metabolite/gene nodes receive messages through HGTConv layers.
    Without reverse edges, only reaction nodes (the dst of all forward edges) get
    updated, leaving metabolite/gene as information dead-ends.
  - GNN+XGBoost separation: avoids multi-objective loss conflict
  - GNN fine-tuning allowed: pretrained encoder must be fine-tunable
  - Edge prediction pretraining: alternative to failed autoencoder
"""

import time
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import HeteroData
from torch_geometric.nn import HGTConv
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import xgboost as xgb


# ---------------------------------------------------------------------------
# 1. Heterogeneous GNN Encoder
# ---------------------------------------------------------------------------

class HGTGNN(nn.Module):
    """
    Heterogeneous Graph Transformer encoder.

    Architecture:
        Input node features → HGTConv(layer 1) → ReLU → Dropout →
        HGTConv(layer 2) → gene-mean-pooling → graph embedding (out_channels)

    Knockout mask is applied by multiplying gene node features:
        gene_feat = [degree, 1] * [1, (1 - knockout)] → knocked-out genes get zeroed
    This preserves the feature dimension while encoding knockout information.
    """

    def __init__(
        self,
        metadata: Tuple[List[str], List[Tuple[str, str, str]]],
        hidden_channels: int = 32,
        out_channels: int = 32,
        num_heads: int = 2,
        num_layers: int = 2,
        dropout: float = 0.1,
        node_feature_dim: int = 2,
    ):
        super().__init__()

        self.node_types = metadata[0]
        self.edge_types = metadata[1]
        self.out_channels = out_channels
        self.num_layers = num_layers
        self.node_feature_dim = node_feature_dim

        # Per-type linear projections to hidden_channels
        self.node_embeddings = nn.ModuleDict()
        for ntype in self.node_types:
            self.node_embeddings[ntype] = nn.Linear(node_feature_dim, hidden_channels)

        # HGTConv layers
        self.convs = nn.ModuleList()
        for i in range(num_layers):
            self.convs.append(
                HGTConv(
                    in_channels=hidden_channels,
                    out_channels=hidden_channels if i < num_layers - 1 else out_channels,
                    metadata=metadata,
                    heads=num_heads,
                )
            )

        self.dropout_layer = nn.Dropout(dropout)

    def forward(
        self,
        x_dict: Dict[str, torch.Tensor],
        edge_index_dict: Dict[Tuple[str, str, str], torch.Tensor],
        knockout_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Forward pass through HGT layers.

        Args:
            x_dict: node features per type
            edge_index_dict: edge indices per type
            knockout_mask: (n_genes,) binary mask, 1=knocked out

        Returns:
            graph_embedding: Tensor of shape (1, out_channels)
        """
        # Project node features to hidden_channels
        h_dict = {}
        h_prev = {}  # Keep previous layer embeddings for None fallback
        for ntype in self.node_types:
            if ntype in x_dict:
                feat = x_dict[ntype]
                # Apply knockout: multiply gene features by (1 - mask)
                if ntype == "gene" and knockout_mask is not None:
                    # Knocked-out genes → features zeroed out
                    scale = (1.0 - knockout_mask).unsqueeze(1)  # (n_genes, 1)
                    feat = feat * scale
                h_dict[ntype] = self.node_embeddings[ntype](feat)

        # HGTConv layers
        for i, conv in enumerate(self.convs):
            h_prev = {k: v.clone() for k, v in h_dict.items()}
            out_dict = conv(h_dict, edge_index_dict)
            # HGTConv returns None for node types that don't receive messages.
            # With reverse edges (reaction→metabolite, reaction→gene),
            # all types receive messages. Fallback to previous embedding
            # for any None outputs (safety net for edge cases).
            h_dict = {}
            for ntype in self.node_types:
                if ntype in out_dict and out_dict[ntype] is not None:
                    h_dict[ntype] = out_dict[ntype]
                elif ntype in h_prev:
                    h_dict[ntype] = h_prev[ntype]
            if i < self.num_layers - 1:
                h_dict = {k: self.dropout_layer(F.relu(v)) for k, v in h_dict.items()}

        # Gene-mean-pooling: average gene node embeddings
        if "gene" in h_dict and h_dict["gene"] is not None:
            graph_emb = h_dict["gene"].mean(dim=0, keepdim=True)
        else:
            embs = [v.mean(dim=0, keepdim=True) for v in h_dict.values() if v is not None]
            graph_emb = torch.cat(embs, dim=1)

        return graph_emb  # (1, out_channels)


# ---------------------------------------------------------------------------
# 2. Edge Prediction Pretraining
# ---------------------------------------------------------------------------

class EdgePredictionPretrainer:
    """
    Pretrain GNN via edge prediction task.

    Randomly masks metabolite-reaction edges and trains the GNN
    to predict whether an edge exists.
    """

    def __init__(
        self,
        model: HGTGNN,
        lr: float = 0.001,
        weight_decay: float = 1e-5,
    ):
        self.model = model
        self.optimizer = torch.optim.Adam(
            model.parameters(), lr=lr, weight_decay=weight_decay
        )

    def train_step(
        self,
        x_dict: Dict[str, torch.Tensor],
        edge_index_dict: Dict[Tuple[str, str, str], torch.Tensor],
        edge_type: Tuple[str, str, str] = None,
        mask_ratio: float = 0.2,
    ) -> float:
        """Single pretraining step: mask edges, predict them."""
        self.model.train()
        self.optimizer.zero_grad()

        # Pick first edge type if not specified
        if edge_type is None:
            edge_type = list(edge_index_dict.keys())[0]

        # Get edge index for the target edge type
        edge_index = edge_index_dict[edge_type]
        n_edges = edge_index.shape[1]

        if n_edges < 2:
            return 0.0

        # Randomly mask some edges
        n_mask = max(1, int(n_edges * mask_ratio))
        perm = torch.randperm(n_edges)
        mask_idx = perm[:n_mask]
        keep_idx = perm[n_mask:]

        # Create masked edge_index_dict
        masked_edge_index_dict = {}
        for etype, eidx in edge_index_dict.items():
            if etype == edge_type:
                masked_edge_index_dict[etype] = eidx[:, keep_idx]
            else:
                masked_edge_index_dict[etype] = eidx

        # Forward pass with masked edges
        h_dict = {}
        for ntype in self.model.node_types:
            if ntype in x_dict:
                h_dict[ntype] = self.model.node_embeddings[ntype](x_dict[ntype])

        for j, conv in enumerate(self.model.convs):
            h_dict = conv(h_dict, masked_edge_index_dict)
            if j < self.model.num_layers - 1:
                h_dict = {k: F.relu(v) for k, v in h_dict.items()}

        # Edge prediction: dot product of source and target node embeddings
        src_type, _, dst_type = edge_type
        if src_type not in h_dict or dst_type not in h_dict:
            return 0.0

        src_embs = h_dict[src_type]
        dst_embs = h_dict[dst_type]

        # Positive edges (masked ones)
        pos_src = edge_index[0, mask_idx]
        pos_dst = edge_index[1, mask_idx]
        pos_scores = (src_embs[pos_src] * dst_embs[pos_dst]).sum(dim=1)

        # Negative edges (random pairs)
        neg_src = torch.randint(0, src_embs.shape[0], (n_mask,))
        neg_dst = torch.randint(0, dst_embs.shape[0], (n_mask,))
        neg_scores = (src_embs[neg_src] * dst_embs[neg_dst]).sum(dim=1)

        # Binary cross-entropy loss
        pos_loss = F.binary_cross_entropy_with_logits(
            pos_scores, torch.ones_like(pos_scores)
        )
        neg_loss = F.binary_cross_entropy_with_logits(
            neg_scores, torch.zeros_like(neg_scores)
        )
        loss = pos_loss + neg_loss

        loss.backward()
        self.optimizer.step()

        return loss.item()


# ---------------------------------------------------------------------------
# 3. GNN + XGBoost Surrogate Pipeline
# ---------------------------------------------------------------------------

class GNNXGBoostSurrogate:
    """
    Full surrogate pipeline: GNN encoder + XGBoost regressor.

    Pipeline:
        1. Pretrain GNN with edge prediction
        2. Extract graph embeddings for all knockout samples
        3. Concatenate embeddings with knockout mask
        4. Train XGBoost on combined features
    """

    def __init__(
        self,
        graph: HeteroData,
        hidden_channels: int = 32,
        out_channels: int = 32,
        num_heads: int = 2,
        num_layers: int = 2,
        dropout: float = 0.1,
        xgb_params: Optional[Dict] = None,
    ):
        self.graph = graph
        self.n_genes = graph["gene"].num_nodes
        self.out_channels = out_channels

        # Extract tensors from graph
        self.x_dict = {
            "metabolite": graph["metabolite"].x,
            "reaction": graph["reaction"].x,
            "gene": graph["gene"].x,
        }
        self.edge_index_dict = {}
        for etype in graph.edge_types:
            self.edge_index_dict[etype] = graph[etype].edge_index

        self.metadata = graph.metadata()

        # Initialize GNN
        self.gnn = HGTGNN(
            metadata=self.metadata,
            hidden_channels=hidden_channels,
            out_channels=out_channels,
            num_heads=num_heads,
            num_layers=num_layers,
            dropout=dropout,
            node_feature_dim=2,
        )

        # Default XGBoost parameters
        self.xgb_params = xgb_params or {
            "objective": "reg:squarederror",
            "max_depth": 6,
            "learning_rate": 0.1,
            "n_estimators": 200,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
            "random_state": 42,
        }

        self.xgb_model = None
        self.is_pretrained = False

    def pretrain(
        self,
        n_epochs: int = 50,
        lr: float = 0.001,
        verbose: bool = True,
    ) -> List[float]:
        """Pretrain GNN with edge prediction task."""
        if verbose:
            print(f"[Pretrain] Edge prediction pretraining for {n_epochs} epochs...")

        pretrainer = EdgePredictionPretrainer(self.gnn, lr=lr)
        losses = []

        for epoch in range(n_epochs):
            # Cycle through all edge types for pretraining
            total_loss = 0.0
            for etype in self.edge_index_dict:
                loss = pretrainer.train_step(
                    self.x_dict, self.edge_index_dict, edge_type=etype
                )
                total_loss += loss
            avg_loss = total_loss / max(len(self.edge_index_dict), 1)
            losses.append(avg_loss)

            if verbose and (epoch + 1) % 10 == 0:
                print(f"  Epoch {epoch+1}/{n_epochs}: avg_loss={avg_loss:.4f}")

        self.is_pretrained = True

        if verbose:
            print(f"  Final loss: {losses[-1]:.4f} (start: {losses[0]:.4f})")
            improvement = (losses[0] - losses[-1]) / max(abs(losses[0]), 1e-8) * 100
            print(f"  Improvement: {improvement:.1f}%")

        return losses

    def _extract_embeddings(
        self,
        knockout_masks: np.ndarray,
    ) -> np.ndarray:
        """
        Extract GNN embeddings for a batch of knockout masks.

        Returns:
            combined_features: (n_samples, out_channels + n_genes)
        """
        self.gnn.eval()
        masks_tensor = torch.tensor(knockout_masks, dtype=torch.float32)
        n_samples = masks_tensor.shape[0]

        all_embeddings = []

        with torch.no_grad():
            for i in range(n_samples):
                ko_mask = masks_tensor[i]  # (n_genes,)
                emb = self.gnn(self.x_dict, self.edge_index_dict, knockout_mask=ko_mask)
                all_embeddings.append(emb.squeeze(0))  # (out_channels,)

        embeddings = torch.stack(all_embeddings).numpy()  # (n_samples, out_channels)

        # Concatenate with knockout mask
        combined = np.concatenate([embeddings, knockout_masks], axis=1)

        return combined  # (n_samples, out_channels + n_genes)

    def fit(
        self,
        knockout_masks: np.ndarray,
        growth_rates: np.ndarray,
        pretrain_epochs: int = 50,
        val_ratio: float = 0.2,
        verbose: bool = True,
    ) -> Dict:
        """
        Train the full surrogate pipeline.
        """
        t0 = time.time()

        # Step 1: Pretrain GNN
        if not self.is_pretrained and pretrain_epochs > 0:
            self.pretrain(n_epochs=pretrain_epochs, verbose=verbose)

        # Step 2: Extract embeddings
        if verbose:
            print(f"[Fit] Extracting GNN embeddings for {len(growth_rates)} samples...")
        t_emb = time.time()
        combined_features = self._extract_embeddings(knockout_masks)
        t_emb = time.time() - t_emb
        if verbose:
            print(f"  Embedding extraction: {t_emb:.1f}s, feature dim={combined_features.shape[1]}")

        # Step 3: Train/test split
        X_train, X_test, y_train, y_test = train_test_split(
            combined_features, growth_rates, test_size=val_ratio, random_state=42
        )

        # Step 4: Train XGBoost
        if verbose:
            print(f"[Fit] Training XGBoost: {X_train.shape[0]} train, {X_test.shape[0]} test")
        self.xgb_model = xgb.XGBRegressor(**self.xgb_params)
        self.xgb_model.fit(X_train, y_train, verbose=False)

        # Step 5: Evaluate
        y_pred_train = self.xgb_model.predict(X_train)
        y_pred_test = self.xgb_model.predict(X_test)

        r2_train = r2_score(y_train, y_pred_train)
        r2_test = r2_score(y_test, y_pred_test)
        rmse_train = np.sqrt(mean_squared_error(y_train, y_pred_train))
        rmse_test = np.sqrt(mean_squared_error(y_test, y_pred_test))
        mae_test = mean_absolute_error(y_test, y_pred_test)

        elapsed = time.time() - t0

        results = {
            "r2_train": r2_train,
            "r2_test": r2_test,
            "rmse_train": rmse_train,
            "rmse_test": rmse_test,
            "mae_test": mae_test,
            "n_train": len(y_train),
            "n_test": len(y_test),
            "feature_dim": combined_features.shape[1],
            "total_time": elapsed,
        }

        if verbose:
            print(f"[Fit] Results:")
            print(f"  R2 train: {r2_train:.4f}, R2 test: {r2_test:.4f}")
            print(f"  RMSE train: {rmse_train:.4f}, RMSE test: {rmse_test:.4f}")
            print(f"  MAE test: {mae_test:.4f}")
            print(f"  Total time: {elapsed:.1f}s")

        return results

    def predict(self, knockout_masks: np.ndarray) -> np.ndarray:
        """Predict growth rates for knockout combinations."""
        if self.xgb_model is None:
            raise RuntimeError("Model not trained. Call fit() first.")
        combined_features = self._extract_embeddings(knockout_masks)
        return self.xgb_model.predict(combined_features)

    def evaluate(
        self,
        knockout_masks: np.ndarray,
        growth_rates: np.ndarray,
        n_folds: int = 5,
    ) -> Dict:
        """Evaluate with k-fold cross-validation."""
        combined_features = self._extract_embeddings(knockout_masks)

        kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
        r2_scores = []
        rmse_scores = []

        for fold, (train_idx, test_idx) in enumerate(kf.split(combined_features)):
            X_train = combined_features[train_idx]
            X_test = combined_features[test_idx]
            y_train = growth_rates[train_idx]
            y_test = growth_rates[test_idx]

            model = xgb.XGBRegressor(**self.xgb_params)
            model.fit(X_train, y_train, verbose=False)

            y_pred = model.predict(X_test)
            r2 = r2_score(y_test, y_pred)
            rmse = np.sqrt(mean_squared_error(y_test, y_pred))

            r2_scores.append(r2)
            rmse_scores.append(rmse)

        results = {
            "r2_mean": np.mean(r2_scores),
            "r2_std": np.std(r2_scores),
            "rmse_mean": np.mean(rmse_scores),
            "rmse_std": np.std(rmse_scores),
            "r2_folds": r2_scores,
            "rmse_folds": rmse_scores,
        }

        print(f"[CV-{n_folds}] R2: {results['r2_mean']:.4f} +/- {results['r2_std']:.4f}")
        print(f"[CV-{n_folds}] RMSE: {results['rmse_mean']:.4f} +/- {results['rmse_std']:.4f}")

        return results


# ---------------------------------------------------------------------------
# 4. XGBoost-only baseline (no GNN)
# ---------------------------------------------------------------------------

class XGBoostOnlyBaseline:
    """
    Baseline: XGBoost with raw knockout mask only (no GNN embeddings).
    This serves as the R2=-0.31 baseline from stages/01.
    """

    def __init__(self, xgb_params: Optional[Dict] = None):
        self.xgb_params = xgb_params or {
            "objective": "reg:squarederror",
            "max_depth": 6,
            "learning_rate": 0.1,
            "n_estimators": 200,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
            "random_state": 42,
        }
        self.xgb_model = None

    def fit(
        self,
        knockout_masks: np.ndarray,
        growth_rates: np.ndarray,
        val_ratio: float = 0.2,
        verbose: bool = True,
    ) -> Dict:
        X_train, X_test, y_train, y_test = train_test_split(
            knockout_masks, growth_rates, test_size=val_ratio, random_state=42
        )

        self.xgb_model = xgb.XGBRegressor(**self.xgb_params)
        self.xgb_model.fit(X_train, y_train, verbose=False)

        y_pred_train = self.xgb_model.predict(X_train)
        y_pred_test = self.xgb_model.predict(X_test)

        results = {
            "r2_train": r2_score(y_train, y_pred_train),
            "r2_test": r2_score(y_test, y_pred_test),
            "rmse_test": np.sqrt(mean_squared_error(y_test, y_pred_test)),
        }

        if verbose:
            print(f"[XGBoost-only] R2 train: {results['r2_train']:.4f}, "
                  f"R2 test: {results['r2_test']:.4f}")

        return results

    def predict(self, knockout_masks: np.ndarray) -> np.ndarray:
        return self.xgb_model.predict(knockout_masks)


# ---------------------------------------------------------------------------
# 5. Test / verification
# ---------------------------------------------------------------------------

def run_tests():
    """Quick verification tests for Module B."""
    print("=" * 60)
    print("Module B: GNN + XGBoost Surrogate -- Verification Tests")
    print("=" * 60)

    all_pass = True

    from module_a_fba_generator import FBAGroundTruthGenerator

    # Test 1: Generate data with Module A
    print("\n--- Test 1: Data Generation (Module A) ---")
    try:
        gen = FBAGroundTruthGenerator()
        data = gen.run(single=True, double=False, random_n=50, n_workers=1)
        print(f"  PASS: Data generated")
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False
        return

    # Test 2: HGTGNN instantiation and forward pass
    print("\n--- Test 2: HGTGNN Forward Pass ---")
    try:
        graph = data["graph"]
        metadata = graph.metadata()

        gnn = HGTGNN(
            metadata=metadata,
            hidden_channels=16,
            out_channels=16,
            num_heads=2,
            num_layers=2,
        )

        x_dict = {
            "metabolite": graph["metabolite"].x,
            "reaction": graph["reaction"].x,
            "gene": graph["gene"].x,
        }
        edge_index_dict = {}
        for etype in graph.edge_types:
            edge_index_dict[etype] = graph[etype].edge_index

        # Standard forward
        emb = gnn(x_dict, edge_index_dict)
        assert emb.shape == (1, 16), f"Embedding shape {emb.shape} != (1, 16)"
        print(f"  Standard forward OK: shape={emb.shape}")

        # Forward with knockout mask
        import torch
        ko_mask = torch.zeros(137, dtype=torch.float32)
        ko_mask[0] = 1.0  # Knock out first gene
        emb_ko = gnn(x_dict, edge_index_dict, knockout_mask=ko_mask)
        assert emb_ko.shape == (1, 16), f"KO embedding shape {emb_ko.shape} != (1, 16)"
        print(f"  Knockout forward OK: shape={emb_ko.shape}")

        # Embeddings should differ when knockout is applied
        diff = (emb - emb_ko).abs().sum().item()
        print(f"  Embedding difference with knockout: {diff:.4f}")
        assert diff > 0, "Knockout should change the embedding"
        print(f"  PASS: Knockout mask affects GNN embedding")
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()
        all_pass = False

    # Test 3: Edge prediction pretraining
    print("\n--- Test 3: Edge Prediction Pretraining ---")
    try:
        surrogate = GNNXGBoostSurrogate(graph, hidden_channels=16, out_channels=16)
        losses = surrogate.pretrain(n_epochs=20, verbose=True)
        assert len(losses) == 20
        print(f"  PASS: Pretraining completed, loss {losses[0]:.4f} -> {losses[-1]:.4f}")
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Test 4: Full surrogate pipeline fit
    print("\n--- Test 4: Full Surrogate Pipeline ---")
    try:
        all_masks = []
        all_growth = []

        for key_mask, key_growth in [
            ("single_ko_mask", "single_ko_growth"),
            ("random_ko_mask", "random_ko_growth"),
        ]:
            if key_mask in data:
                mask = data[key_mask]
                if isinstance(mask, torch.Tensor):
                    mask = mask.numpy()
                all_masks.append(mask)
                all_growth.append(data[key_growth])

        X = np.concatenate(all_masks, axis=0)
        y = np.concatenate(all_growth, axis=0)

        print(f"  Dataset: {X.shape[0]} samples, {X.shape[1]} genes")
        print(f"  Growth range: [{y.min():.4f}, {y.max():.4f}], mean={y.mean():.4f}")

        surrogate = GNNXGBoostSurrogate(graph, hidden_channels=16, out_channels=16)
        results = surrogate.fit(X, y, pretrain_epochs=10, verbose=True)

        print(f"  R2 test: {results['r2_test']:.4f}")
        print(f"  PASS: Full pipeline completed")
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()
        all_pass = False

    # Test 5: XGBoost-only baseline
    print("\n--- Test 5: XGBoost-only Baseline ---")
    try:
        baseline = XGBoostOnlyBaseline()
        baseline_results = baseline.fit(X, y, verbose=True)
        print(f"  PASS: Baseline R2 test = {baseline_results['r2_test']:.4f}")
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Test 6: Cross-validation
    print("\n--- Test 6: 3-Fold Cross-Validation ---")
    try:
        cv_results = surrogate.evaluate(X, y, n_folds=3)
        print(f"  PASS: CV R2 = {cv_results['r2_mean']:.4f} +/- {cv_results['r2_std']:.4f}")
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Summary
    print("\n" + "=" * 60)
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED -- see above")
    print("=" * 60)


if __name__ == "__main__":
    run_tests()
