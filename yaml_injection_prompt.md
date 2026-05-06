# Prompt: Inject Research Object YAML Layer Into Existing Claude-Code Research Agent

You are working inside an existing research-agent project that already has:
- CLAUDE.md
- docs/01_literature_review.md
- docs/02_framing.md
- docs/03_planning.md
- docs/04_analysis.md
- docs/05_interpretation.md
- docs/06_git_policy.md
- docs/07_experiment_failure_reports.md
- stages/
- outputs/

The user has provided a ZIP bundle containing YAML templates under:

objects/current/*.yaml

Your task is to add these YAML research objects to the project and integrate them into the existing workflow conceptually, without rewriting the current stage guides.

## Do not modify

Do not edit the following files unless the user explicitly asks later:
- CLAUDE.md
- docs/01_literature_review.md
- docs/02_framing.md
- docs/03_planning.md
- docs/04_analysis.md
- docs/05_interpretation.md
- docs/06_git_policy.md
- docs/07_experiment_failure_reports.md
- stages/*
- outputs/*

## Goal

Add a lightweight human-auditable research object layer.

The existing workflow tells the agent what to do at each research stage.
The new YAML objects should track the current decision state of the research project:
- idea abstraction
- novelty status
- experiment contract
- metric validity
- run result summary
- validation readiness
- pivot diagnosis
- allowed claim strength

These objects do not replace stages/.
They complement stages/.

Use this distinction:

- stages/ = refined narrative state
- outputs/ = detailed run artifacts, logs, reports, code, results
- objects/current/ = compact current decision state for human review
- objects/history/ = archived snapshots of object states before major changes

## Step 1. Inspect project structure

Read:
1. CLAUDE.md
2. docs/01_literature_review.md
3. docs/02_framing.md
4. docs/03_planning.md
5. docs/04_analysis.md
6. docs/05_interpretation.md
7. docs/06_git_policy.md
8. docs/07_experiment_failure_reports.md

Then report:
- current stage
- existing stage flow
- where objects/current/ will sit
- whether objects/current/ already exists

Do not perform experiments.

## Step 2. Install YAML templates

Create directories if missing:

```text
objects/current/
objects/history/
```

Copy the YAML templates from the provided ZIP into:

```text
objects/current/
```

Expected files:

```text
objects/current/idea_abstraction_card.yaml
objects/current/novelty_ledger.yaml
objects/current/experiment_contract.yaml
objects/current/evaluation_validity_card.yaml
objects/current/result_card.yaml
objects/current/validation_readiness_card.yaml
objects/current/pivot_diagnosis_card.yaml
objects/current/claim_card.yaml
```

Also create:

```text
objects/history/.gitkeep
```

Do not overwrite existing object files unless they are empty templates.
If an object file already exists and contains project-specific content, create a backup under:

```text
objects/history/pre_injection_backup_YYYYMMDD_HHMMSS/
```

## Step 3. Validate YAML syntax

For each YAML file:
- check that it parses as valid YAML
- check that it has object_type
- check that it has status
- check that it has linked_stage
- check that it has evidence
- check that it has uncertainty
- check that it has human_review_needed

If validation fails:
- fix only YAML syntax
- do not fill in project-specific facts unless directly supported by existing stages/ or outputs/

## Step 4. Create a lightweight README

Create:

```text
objects/README.md
```

The README must explain:

1. Purpose of objects/current/
2. Difference between stages/, outputs/, and objects/
3. When each object should be updated
4. Required object checks before Analysis
5. Required object checks after Analysis
6. Pivot rule: PivotDiagnosisCard is required before changing research direction

Keep this README concise.

## Step 5. Do not pre-fill unsupported claims

Do not invent project-specific results.
Use existing stages/ and outputs/ only if they contain explicit evidence.

If uncertain, leave fields as:
- "TBD"
- null
- []
- "unknown"

## Step 6. Report completion

Return a summary:

### Added files
List all files created.

### Backups
State whether any existing object files were backed up.

### YAML validation
State whether all object files parse correctly.

### Integration notes
Briefly explain how to use the YAML layer:
- before a run
- before Analysis
- after a run
- before pivot
- before writing claims

### Recommended next prompt
Suggest that the next step is to run a novelty audit or a pre-analysis gate, depending on current stage.
