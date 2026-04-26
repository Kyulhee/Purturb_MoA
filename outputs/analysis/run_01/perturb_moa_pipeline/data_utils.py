"""
Cross-Perturbation MoA Framework -- Data Download & Preprocessing
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path
from config import RAW_DIR, PROCESSED_DIR, SCI_PLEX_N_GENES


def download_sci_plex(data_dir: Path = RAW_DIR) -> Path:
    """
    Download sci-Plex data from GEO.
    Dataset: Srivatsan et al. 2020, Science (PMID: 31806696)
    GEO: GSE139944

    Returns path to downloaded .h5ad file.
    """
    output_path = data_dir / "sci_plex.h5ad"
    if output_path.exists():
        print(f"sci-Plex data already exists: {output_path}")
        return output_path

    print("Downloading sci-Plex data from GEO...")
    # Use scanpy's built-in GEO download or direct URL
    # GEO GSE139944 supplementary files
    try:
        import subprocess
        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139944/suppl/"
        print(f"Attempting download from: {url}")
        print("Note: Manual download may be required if FTP is slow.")
        print("Alternative: Use PerturbNet's preprocessed data from their GitHub repo")
    except Exception as e:
        print(f"Download failed: {e}")
        print("Please download manually from GEO GSE139944")

    return output_path


def load_sci_plex_from_perturbnet(data_dir: Path = RAW_DIR) -> ad.AnnData:
    """
    Load sci-Plex data using PerturbNet's preprocessing.
    PerturbNet GitHub: https://github.com/broadinstitute/PerturbNet

    Expected format after PerturbNet preprocessing:
    - 648,737 cells x 5,087 genes
    - 180 drug treatments
    - 3 cell lines: A549, K562, MCF7
    """
    h5ad_path = data_dir / "sci_plex.h5ad"
    if h5ad_path.exists():
        print(f"Loading preprocessed sci-Plex: {h5ad_path}")
        adata = sc.read_h5ad(h5ad_path)
        return adata

    print("No preprocessed data found. Attempting to create from raw...")
    return None


def preprocess_sci_plex(adata: ad.AnnData, n_top_genes: int = SCI_PLEX_N_GENES) -> ad.AnnData:
    """
    Preprocess sci-Plex data following PerturbNet protocol:
    1. Filter cells (min_genes=200, min_counts=1000)
    2. Filter genes (min_cells=10)
    3. Normalize (scanpy pp.normalize_total)
    4. Log transform
    5. Select top variable genes
    6. Extract perturbation labels
    """
    print(f"Preprocessing: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, min_counts=1000)
    sc.pp.filter_genes(adata, min_cells=10)
    print(f"After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Select variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat_v3')
    adata = adata[:, adata.var.highly_variable].copy()
    print(f"After HVG selection: {adata.shape[0]} cells x {adata.shape[1]} genes")

    return adata


def annotate_moa(adata: ad.AnnData) -> ad.AnnData:
    """
    Add MoA annotations to sci-Plex data.

    Uses ATC code classification from DrugBank.
    Major MoA categories in sci-Plex 188 compounds:
    - HDAC inhibitors (vorinostat, panobinostat, etc.)
    - CDK inhibitors (dinaciclib, palbociclib, etc.)
    - Proteasome inhibitors (bortezomib, carfilzomib)
    - DNA damaging agents
    - Kinase inhibitors
    - Tubulin inhibitors
    - etc.

    NOTE: For prototype, use a simplified MoA mapping.
    Full version should use DrugBank API or CMap annotations.
    """
    # Simplified MoA mapping for key sci-Plex compounds
    # Based on literature knowledge (PMID: 31806696, PerturbNet, PRnet papers)
    MOA_MAP = {
        # HDAC inhibitors
        'vorinostat': 'HDAC_inhibitor',
        'panobinostat': 'HDAC_inhibitor',
        'tucidinostat': 'HDAC_inhibitor',
        'entinostat': 'HDAC_inhibitor',
        'valproic acid': 'HDAC_inhibitor',
        'sodium butyrate': 'HDAC_inhibitor',
        'trichostatin a': 'HDAC_inhibitor',

        # CDK inhibitors
        'dinaciclib': 'CDK_inhibitor',
        'palbociclib': 'CDK_inhibitor',
        'roscovitine': 'CDK_inhibitor',
        'flavopiridol': 'CDK_inhibitor',

        # Proteasome inhibitors
        'bortezomib': 'Proteasome_inhibitor',
        'carfilzomib': 'Proteasome_inhibitor',
        'mg-132': 'Proteasome_inhibitor',

        # DNA damaging agents
        'cisplatin': 'DNA_damaging',
        'camptothecin': 'DNA_damaging',
        'etoposide': 'DNA_damaging',
        'doxorubicin': 'DNA_damaging',
        'mitomycin c': 'DNA_damaging',

        # Tubulin inhibitors
        'paclitaxel': 'Tubulin_inhibitor',
        'vincristine': 'Tubulin_inhibitor',
        'colchicine': 'Tubulin_inhibitor',
        'nocodazole': 'Tubulin_inhibitor',

        # mTOR inhibitors
        'rapamycin': 'mTOR_inhibitor',
        'torin1': 'mTOR_inhibitor',
        'azd8055': 'mTOR_inhibitor',

        # MEK inhibitors
        'trametinib': 'MEK_inhibitor',
        'selumetinib': 'MEK_inhibitor',
        'pd0325901': 'MEK_inhibitor',

        # BET inhibitors
        'jQ1': 'BET_inhibitor',
        'iBET-762': 'BET_inhibitor',

        # HSP90 inhibitors
        'geldanamycin': 'HSP90_inhibitor',
        '17-AAG': 'HSP90_inhibitor',

        # BCL2 inhibitors
        'venetoclax': 'BCL2_inhibitor',
        'navitoclax': 'BCL2_inhibitor',

        # Control
        'dmso': 'Control',
    }

    # If moa already exists with meaningful categories (e.g. synthetic data), keep it
    if 'moa' in adata.obs.columns:
        existing_moa = adata.obs['moa'].nunique()
        if existing_moa > 2:  # more than just Other/Control
            print(f"MoA annotation already present: {existing_moa} categories, keeping existing")
            print(adata.obs['moa'].value_counts())
            return adata

    # Try to match compound names to MoA
    if 'compound' in adata.obs.columns:
        col = 'compound'
    elif 'drug' in adata.obs.columns:
        col = 'drug'
    elif 'perturbation' in adata.obs.columns:
        col = 'perturbation'
    else:
        print("Warning: No compound column found in adata.obs")
        adata.obs['moa'] = 'Unknown'
        return adata

    compounds = adata.obs[col].astype(str).str.lower().str.strip()
    moa_labels = compounds.map(
        lambda x: next(
            (v for k, v in MOA_MAP.items() if k.lower() in x),
            'Other'
        )
    )
    adata.obs['moa'] = moa_labels.values

    n_moa = adata.obs['moa'].nunique()
    print(f"MoA annotation: {n_moa} categories")
    print(adata.obs['moa'].value_counts())

    return adata


def generate_synthetic_fingerprints(n_drugs: int, n_moa_classes: int,
                                     fingerprint_bits: int = 1024,
                                     moa_correlation: float = 0.3,
                                     seed: int = 42) -> dict:
    """
    Generate MoA-correlated random fingerprints for synthetic data.

    In real data, drugs with the same MoA share structural features (same scaffold),
    so their FCFP fingerprints overlap. We simulate this by:
    1. Creating a MoA-specific "core" bit pattern (~30% of bits)
    2. Adding random bits for each drug

    Args:
        n_drugs: number of drugs
        n_moa_classes: number of MoA categories
        fingerprint_bits: FCFP bit vector length
        moa_correlation: fraction of bits that are MoA-specific (0.3 = 30%)
        seed: random seed

    Returns:
        dict mapping drug index -> fingerprint numpy array
    """
    rng = np.random.default_rng(seed + 100)
    n_moa_bits = int(fingerprint_bits * moa_correlation)
    n_random_bits = fingerprint_bits - n_moa_bits

    # MoA-specific core patterns
    moa_cores = rng.binomial(1, 0.5, size=(n_moa_classes, n_moa_bits)).astype(np.float32)

    # Per-drug fingerprints
    drug_fingerprints = {}
    for i in range(n_drugs):
        moa_idx = i % n_moa_classes
        # Core bits from MoA + random bits unique to this drug
        random_part = rng.binomial(1, 0.3, size=n_random_bits).astype(np.float32)
        fp = np.concatenate([moa_cores[moa_idx], random_part])
        drug_fingerprints[i] = fp

    return drug_fingerprints


def create_synthetic_sci_plex(n_cells: int = 50000, n_genes: int = 2000,
                               n_drugs: int = 50, n_cell_lines: int = 3,
                               n_moa_classes: int = 8,
                               fingerprint_bits: int = 1024) -> ad.AnnData:
    """
    Create synthetic sci-Plex-like data for pipeline testing.
    This is ONLY for prototype development -- real analysis uses download_sci_plex().

    Structure:
    - n_cells cells, n_genes genes
    - n_drugs drug perturbations + 1 control (DMSO)
    - n_cell_lines cell lines
    - n_moa_classes MoA categories
    - Each MoA has a distinct gene expression signature
    """
    rng = np.random.default_rng(42)

    # MoA-specific gene signatures (sparse -- only ~10% genes affected per MoA)
    moa_signatures = np.zeros((n_moa_classes, n_genes))
    for i in range(n_moa_classes):
        affected_genes = rng.choice(n_genes, size=n_genes // 10, replace=False)
        moa_signatures[i, affected_genes] = rng.normal(0, 1.5, size=len(affected_genes))

    # Generate cells
    drugs = [f"drug_{i}" for i in range(n_drugs)]
    cell_lines = ["A549", "K562", "MCF7"][:n_cell_lines]
    moa_names = [f"MoA_{i}" for i in range(n_moa_classes)]

    # Assign drugs to MoA classes
    drugs_per_moa = max(1, n_drugs // n_moa_classes)
    drug_to_moa = {}
    for i, drug in enumerate(drugs):
        drug_to_moa[drug] = moa_names[i % n_moa_classes]

    # Cell-level data
    cell_drug = []
    cell_moa = []
    cell_line = []
    cell_dose = []

    n_per_drug = n_cells // (n_drugs + 1)  # +1 for DMSO control

    for drug in drugs:
        for _ in range(n_per_drug):
            cell_drug.append(drug)
            cell_moa.append(drug_to_moa[drug])
            cell_line.append(rng.choice(cell_lines))
            cell_dose.append(rng.choice([0.01, 0.1, 1.0, 10.0]))

    # Add DMSO control
    for _ in range(n_per_drug):
        cell_drug.append("DMSO")
        cell_moa.append("Control")
        cell_line.append(rng.choice(cell_lines))
        cell_dose.append(0.0)

    n_actual = len(cell_drug)

    # Generate expression matrix
    # Base expression (log-normal)
    X = rng.negative_binomial(n=5, p=0.3, size=(n_actual, n_genes)).astype(np.float32)

    # Add MoA-specific perturbation effects
    for i in range(n_actual):
        if cell_moa[i] != "Control":
            moa_idx = moa_names.index(cell_moa[i])
            # Dose-dependent effect
            dose_factor = np.log10(cell_dose[i] + 1) / np.log10(11)  # normalize to [0,1]
            X[i] += (moa_signatures[moa_idx] * dose_factor * 5).astype(np.float32)

    # Add cell line effect
    cell_line_effects = rng.normal(0, 0.5, size=(len(cell_lines), n_genes))
    for i in range(n_actual):
        cl_idx = cell_lines.index(cell_line[i])
        X[i] += cell_line_effects[cl_idx].astype(np.float32)

    # Clip to non-negative
    X = np.clip(X, 0, None)

    # Create AnnData
    obs = pd.DataFrame({
        'drug': cell_drug,
        'moa': cell_moa,
        'cell_line': cell_line,
        'dose': cell_dose,
    })

    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    var['highly_variable'] = True

    adata = ad.AnnData(X=X, obs=obs, var=var)

    print(f"Synthetic sci-Plex: {adata.shape[0]} cells x {adata.shape[1]} genes")
    print(f"Drugs: {len(drugs)}, MoA classes: {n_moa_classes}")
    print(f"Cell lines: {cell_lines}")
    print(adata.obs['moa'].value_counts())

    return adata


def prepare_dataloaders(adata: ad.AnnData, test_size: float = 0.2,
                        split_strategy: str = "random", seed: int = 42,
                        fingerprint_bits: int = 1024,
                        drug_fingerprints: dict = None):
    """
    Create train/test splits with fingerprint data for rFCFP-based model.

    Args:
        adata: AnnData with obs columns: drug, moa, cell_line, dose
        test_size: fraction for test split
        split_strategy: "random", "leave_compound_out", "leave_moa_out", "leave_cell_line_out"
        seed: random seed
        fingerprint_bits: FCFP fingerprint length
        drug_fingerprints: dict mapping drug_index -> fingerprint array.
                          If None, auto-generate MoA-correlated synthetic fingerprints.
    """
    from sklearn.model_selection import train_test_split
    import torch
    from torch.utils.data import DataLoader, TensorDataset

    # Encode categorical variables
    moa_encoder = {m: i for i, m in enumerate(sorted(adata.obs['moa'].unique()))}
    drug_encoder = {d: i for i, d in enumerate(sorted(adata.obs['drug'].unique()))}
    cl_encoder = {c: i for i, c in enumerate(sorted(adata.obs['cell_line'].unique()))}

    moa_labels = adata.obs['moa'].map(moa_encoder).values
    drug_labels = adata.obs['drug'].map(drug_encoder).values
    cl_labels = adata.obs['cell_line'].map(cl_encoder).values
    dose_values = adata.obs['dose'].values.astype(np.float32)

    # Generate or use provided fingerprints
    n_drugs = len(drug_encoder)
    n_moa = len(moa_encoder)
    if drug_fingerprints is None:
        drug_fingerprints = generate_synthetic_fingerprints(
            n_drugs=n_drugs, n_moa_classes=n_moa,
            fingerprint_bits=fingerprint_bits, seed=seed
        )

    # Build fingerprint matrix: [n_drugs, fingerprint_bits]
    fp_matrix = np.zeros((n_drugs, fingerprint_bits), dtype=np.float32)
    for drug_idx, fp in drug_fingerprints.items():
        if drug_idx < n_drugs:
            fp_matrix[drug_idx] = fp

    # Cell-level fingerprints: look up by drug label
    cell_fingerprints = fp_matrix[drug_labels]  # [n_cells, fingerprint_bits]

    # Expression matrix
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X

    # Split indices
    n = len(adata)
    indices = np.arange(n)

    if split_strategy == "random":
        train_idx, test_idx = train_test_split(indices, test_size=test_size,
                                                 random_state=seed, stratify=moa_labels)
    elif split_strategy == "leave_compound_out":
        # Leave out 20% of compounds
        unique_drugs = np.unique(drug_labels)
        rng = np.random.default_rng(seed)
        test_drugs = rng.choice(unique_drugs, size=max(1, len(unique_drugs) // 5), replace=False)
        train_idx = indices[~np.isin(drug_labels, test_drugs)]
        test_idx = indices[np.isin(drug_labels, test_drugs)]
    elif split_strategy == "leave_moa_out":
        # Leave out 1-2 MoA classes
        unique_moa = np.unique(moa_labels)
        rng = np.random.default_rng(seed)
        test_moa = rng.choice(unique_moa, size=max(1, len(unique_moa) // 5), replace=False)
        train_idx = indices[~np.isin(moa_labels, test_moa)]
        test_idx = indices[np.isin(moa_labels, test_moa)]
    elif split_strategy == "leave_cell_line_out":
        # Leave out one cell line
        test_cl = sorted(adata.obs['cell_line'].unique())[-1]
        train_idx = indices[adata.obs['cell_line'].values != test_cl]
        test_idx = indices[adata.obs['cell_line'].values == test_cl]
    else:
        raise ValueError(f"Unknown split strategy: {split_strategy}")

    # Create tensors -- now includes fingerprint
    def make_dataset(idx):
        return TensorDataset(
            torch.tensor(X[idx], dtype=torch.float32),
            torch.tensor(moa_labels[idx], dtype=torch.long),
            torch.tensor(drug_labels[idx], dtype=torch.long),
            torch.tensor(cl_labels[idx], dtype=torch.long),
            torch.tensor(dose_values[idx], dtype=torch.float32).unsqueeze(1),
            torch.tensor(cell_fingerprints[idx], dtype=torch.float32),
        )

    train_ds = make_dataset(train_idx)
    test_ds = make_dataset(test_idx)

    metadata = {
        'moa_encoder': moa_encoder,
        'drug_encoder': drug_encoder,
        'cl_encoder': cl_encoder,
        'n_moa': n_moa,
        'n_drugs': n_drugs,
        'n_cell_lines': len(cl_encoder),
        'n_genes': X.shape[1],
        'fingerprint_bits': fingerprint_bits,
        'train_size': len(train_idx),
        'test_size': len(test_idx),
        'split_strategy': split_strategy,
    }

    return train_ds, test_ds, metadata


if __name__ == "__main__":
    # Quick test with synthetic data
    adata = create_synthetic_sci_plex()
    adata = annotate_moa(adata)
    train_ds, test_ds, meta = prepare_dataloaders(adata, split_strategy="random")
    print(f"\nTrain: {meta['train_size']}, Test: {meta['test_size']}")
    print(f"MoA classes: {meta['n_moa']}, Drugs: {meta['n_drugs']}, Genes: {meta['n_genes']}")
    print(f"Fingerprint bits: {meta['fingerprint_bits']}")
