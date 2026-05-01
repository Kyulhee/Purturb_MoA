# BioEval Run 20 Report — GEARS Deep Learning Model (B2)

**Date**: 2026-05-02
**Runtime**: 577.5s
**Script**: `run_20_gears_training.py`
**Environment**: ai_env (Python 3.11, PyTorch 2.5.1, PyG 2.7.0, GEARS)

## Objective

Obtain real DL model (GEARS) predictions for H3 'DL > baseline under BioEval' testing.

## Norman

| Model | R2 | Dir_deg | DEG_auprc | f1@50 | MSE |
|-------|----|---------|-----------|-------|-----|
| GEARS | FAILED | | | | |
| Ridge | 0.929 | 1.000 | 0.949 | 0.560 | 0.0002 |
| mean_predictor | -0.004 | 0.000 | 0.022 | 0.017 | 0.0053 |
| mean_effect | -0.123 | 0.833 | 0.423 | 0.334 | 0.0040 |

## K562

| Model | R2 | Dir_deg | DEG_auprc | f1@50 | MSE |
|-------|----|---------|-----------|-------|-----|
| GEARS | FAILED | | | | |
| Ridge | 0.611 | 0.989 | 0.634 | 0.381 | 0.0022 |
| mean_predictor | -0.027 | 0.000 | 0.036 | 0.022 | 0.0082 |
| mean_effect | 0.038 | 0.789 | 0.186 | 0.150 | 0.0070 |

## RPE1

| Model | R2 | Dir_deg | DEG_auprc | f1@50 | MSE |
|-------|----|---------|-----------|-------|-----|
| GEARS | FAILED | | | | |
| Ridge | 0.695 | 0.983 | 0.725 | 0.281 | 0.0036 |
| mean_predictor | -0.027 | 0.000 | 0.071 | 0.011 | 0.0201 |
| mean_effect | -0.148 | 0.874 | 0.517 | 0.205 | 0.0118 |
