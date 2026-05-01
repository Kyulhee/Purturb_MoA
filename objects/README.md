# objects/ — Research Object Layer

## Purpose

`objects/current/` tracks the **compact current decision state** of the research project. It answers "where are we now?" at a glance, without reading through full stage narratives or run logs.

## Three-Layer Distinction

| Layer | Location | What it holds | Update frequency |
|-------|----------|---------------|-----------------|
| **stages/** | Project root | Refined narrative state — design decisions and their rationale | At stage transitions |
| **outputs/** | Project root | Detailed run artifacts — code, logs, results, reports | Every run |
| **objects/** | This directory | Compact decision state — current status of key research dimensions | Before/after runs and pivots |

## Objects Overview

| Object | Linked stage | Updated when |
|--------|-------------|--------------|
| `idea_abstraction_card.yaml` | literature_review | Literature search completes or gap reassessed |
| `novelty_ledger.yaml` | framing | Research question, novelty assessment, or competitor landscape changes |
| `experiment_contract.yaml` | planning | Experiment design finalized; must be confirmed before analysis |
| `evaluation_validity_card.yaml` | planning | Metrics, baselines, or validation strategy defined |
| `result_card.yaml` | analysis | Run completes and results are interpreted |
| `validation_readiness_card.yaml` | analysis | Evidence status assessed against claim requirements |
| `pivot_diagnosis_card.yaml` | analysis | **Required before any research direction change** |
| `claim_card.yaml` | interpretation | Claims are drafted, strengthened, or weakened |

## Required Checks

### Before Analysis
- [ ] `experiment_contract.yaml` — hypothesis and success criteria defined
- [ ] `evaluation_validity_card.yaml` — metric-construct alignment checked
- [ ] `novelty_ledger.yaml` — novelty status still holds

### After Analysis (each run)
- [ ] `result_card.yaml` — results filled with actual metrics
- [ ] `validation_readiness_card.yaml` — evidence status updated

### Before Pivot
- [ ] `pivot_diagnosis_card.yaml` — **mandatory** before any direction change or loopback

### Before Writing Claims
- [ ] `claim_card.yaml` — evidence mapped, attack points identified
- [ ] `validation_readiness_card.yaml` — claim strength justified by evidence

## History

`objects/history/` stores archived snapshots of object states before major changes. When an object with non-trivial content is about to be overwritten, archive it first:

```
objects/history/pre_injection_backup_YYYYMMDD_HHMMSS/
```

## Rules

1. Never overwrite an object with existing content without archiving to `history/`
2. Never fill unsupported claims — use `TBD`, `null`, `[]`, or `unknown`
3. Each object must have: `object_type`, `status`, `linked_stage`, `evidence`, `uncertainty`, `human_review_needed`
