# BioEval Run 20b Report — Norman GEARS (GO graph fix)

**Date**: 2026-05-02
**Runtime**: 1514.5s

## Norman GEARS Results

| Model | R2 | Dir_deg | DEG_auprc | f1@50 | MSE |
|-------|----|---------|-----------|-------|-----|
| GEARS | -0.699 | 0.422 | 0.183 | 0.156 | 0.0074 |
| Ridge | 0.896 | 1.000 | 0.925 | 0.551 | 0.0004 |
| mean_predictor | -0.004 | 0.000 | 0.022 | 0.017 | 0.0053 |

### H3: GEARS vs Baselines

- vs Ridge: GEARS wins 0/4 metrics
- vs mean_predictor: GEARS wins 3/4 metrics