"""
Cross-Perturbation MoA Framework -- Configuration
"""
import os
from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
MODEL_DIR = PROJECT_ROOT / "models"

# Create directories
for d in [RAW_DIR, PROCESSED_DIR, RESULTS_DIR, MODEL_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Data sources
SCI_PLEX_GEO = "GSE139944"  # sci-Plex dataset
NORMAN_GEO = "GSE133344"    # Norman 2019 (approximate, verify)
REPLOGLE_GEO = None          # To be confirmed

# sci-Plex parameters (from PerturbNet paper)
SCI_PLEX_N_CELLS = 648737
SCI_PLEX_N_GENES = 5087
SCI_PLEX_N_DRUGS = 180       # after PerturbNet filtering
SCI_PLEX_CELL_LINES = ["A549", "K562", "MCF7"]

# Model hyperparameters
class ModelConfig:
    # Latent space
    latent_dim = 128
    basal_dim = 64
    perturbation_dim = 32
    covariate_dim = 16

    # GenotypeVAE
    genotype_hidden_dims = [512, 256]
    genotype_latent_dim = 128

    # Drug Encoder (rFCFP)
    fingerprint_bits = 1024     # FCFP fingerprint length
    drug_hidden_dims = [256, 128]  # MLP hidden dims for fingerprint encoder
    drug_latent_dim = 32       # Must match perturbation_dim

    # Metric learning / zero-shot MoA
    use_metric_learning = True  # Enable prototype-based zero-shot MoA
    prototype_update_momentum = 0.9  # EMA for class prototypes

    # Decoder
    decoder_hidden_dims = [256, 512]

    # Contrastive loss
    temperature = 0.1
    lambda_recon = 1.0
    lambda_adversarial = 0.5
    lambda_contrastive = 0.3
    lambda_classification = 0.5

    # Training
    batch_size = 256
    learning_rate = 1e-4
    n_epochs = 100
    early_stop_patience = 10

    # Evaluation
    eval_every = 5
    top_k = [1, 3, 5]

# Random seed
SEED = 42
