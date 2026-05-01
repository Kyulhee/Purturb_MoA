"""Debug: check K562 data structure for gene-level feature design"""
import scanpy, numpy as np, os, sys

DATA_DIR = r"C:\test\llm_project\nexus-science-win\outputs\analysis\run_04\data"
path = os.path.join(DATA_DIR, "gears_data", "replogle_k562_essential", "perturb_processed.h5ad")
adata = scanpy.read_h5ad(path)
ctrl_mask = adata.obs["condition"] == "ctrl"
ctrl_data = adata[ctrl_mask].X
if hasattr(ctrl_data, "toarray"):
    ctrl_data = ctrl_data.toarray()
ctrl_data = np.asarray(ctrl_data, dtype=np.float32)
print(f"shape={ctrl_data.shape}, dtype={ctrl_data.dtype}")
print(f"NaN={np.isnan(ctrl_data).sum()}, Inf={np.isinf(ctrl_data).sum()}")
print(f"min={ctrl_data.min():.4f}, max={ctrl_data.max():.4f}, mean={ctrl_data.mean():.4f}")

# Try PCA with just 10 components
from sklearn.decomposition import PCA
pca = PCA(n_components=10, random_state=42)
ctrl_pca = pca.fit_transform(ctrl_data)
print(f"PCA OK. Explained var: {pca.explained_variance_ratio_.round(4)}")

# Project a perturbation
pert = "AAMP+ctrl"
pert_mask = adata.obs["condition"] == pert
pert_data = adata[pert_mask].X
if hasattr(pert_data, "toarray"):
    pert_data = pert_data.toarray()
pert_data = np.asarray(pert_data, dtype=np.float32)
pert_pca = pca.transform(pert_data)
shift = pert_pca.mean(axis=0) - ctrl_pca.mean(axis=0)
print(f"Pert {pert}: {pert_mask.sum()} cells, PCA shift={shift[:5].round(3)}")
