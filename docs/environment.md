# Environment Reference

## Two Python Environments

### 1. System Python (analysis runtime)
- **Path**: `C:\Users\hgh97\AppData\Local\Python\pythoncore-3.14-64\`
- **Version**: Python 3.14.4
- **CUDA**: torch 2.11.0+cu128, CUDA=True (RTX 4060 Ti 8GB)
- **Usage**: BioEval analysis scripts (run_13–run_19), data processing, statistical tests
- **Key packages**: numpy 2.4.4, scipy, sklearn, scanpy, pandas, cell-gears (torch-based)

### 2. ai_env (deep learning)
- **Path**: `C:\Users\hgh97\miniconda3\envs\ai_env\`
- **Version**: Python 3.11.15
- **CUDA**: torch 2.5.1, CUDA=True (RTX 4060 Ti 8GB)
- **Usage**: GEARS/CPA training, PyG-dependent models
- **Key packages**: torch 2.5.1, PyG 2.7.0, scanpy 1.11.5, anndata 0.12.10, sklearn 1.8.0, scipy 1.17.1, numpy 2.0.1

### When to use which
- **System Python**: `python script.py` — analysis scripts, Ridge LOO, bootstrap, metrics
- **ai_env**: `"C:/Users/hgh97/miniconda3/envs/ai_env/python.exe" script.py` — GEARS training, GNN models

### Hardware
- GPU: NVIDIA GeForce RTX 4060 Ti, 8188 MiB VRAM
- Driver: CUDA 12.x (cu128 compatible)

### Data Locations
- GEARS benchmark data: `outputs/analysis/run_04/data/gears_data/`
  - `norman/` — Norman 2019 h5ad
  - `replogle_k562_essential/` — Replogle K562 h5ad
  - `replogle_rpe1_essential/` — Replogle RPE1 h5ad
  - `essential_all_data_pert_genes.pkl` — perturbation gene lists
  - `gene2go_all.pkl` — GO pathway annotations
