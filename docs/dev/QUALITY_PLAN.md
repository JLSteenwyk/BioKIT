# BioKIT Quality Plan (2 Weeks)

## Goal
Improve maintainability, reliability, and release safety without destabilizing CLI behavior.

## Principles
- Keep changes PR-sized and reversible.
- Preserve all current command outputs unless explicitly versioned.
- Add quality gates only after the first passing baseline.

## Week 1

### PR 1: Tooling Baseline (completed in progress)
- Add `requirements-dev.txt` for local/CI developer tooling.
- Add `pyproject.toml` tool configuration (`mypy`, `ruff`).
- Add `.pre-commit-config.yaml`.
- Add Make targets: `lint`, `typecheck`, `quality`, `precommit.*`.
- Add CI `quality` job (`lint` + `typecheck`).

Acceptance:
- `make lint` passes.
- `make typecheck` passes for `biokit/helpers`.
- Existing unit/integration tests pass.

### PR 2: Type-Safety Pilot Expansion
- Expand `mypy` scope from `biokit/helpers` to:
  - `biokit/services/coding_sequences`
  - `biokit/services/alignment`
- Add focused annotations and remove obvious `Any` leaks.
- Keep `ignore_missing_imports = true` while baseline matures.

Acceptance:
- `make typecheck` passes on expanded modules.
- No regression in integration test suite.

### PR 3: Private API Removal (high-churn files first)
- Replace `_seq` / `_data` use with public Biopython APIs in:
  - `biokit/services/coding_sequences/gene_wise_relative_synonymous_codon_usage.py`
  - `biokit/helpers/calculate_character_frequencies.py` (already done)
  - `biokit/services/text/*` where direct internal fields are used

Acceptance:
- Same output for current integration fixtures.
- No new deprecation warnings from touched code paths.

## Week 2

### PR 4: Packaging Modernization
- Migrate metadata from `setup.py` to `pyproject.toml` (`[project]`).
- Keep `setup.py` shim temporarily for compatibility, then remove.
- Separate runtime vs dev/documentation dependencies.

Acceptance:
- `pip install .` and `pip install -e .` both work.
- CLI entry points match current behavior.

### PR 5: Test and CI Hardening
- Add a matrix for supported Python versions in CI quality + tests.
- Add regression tests for CLI exit codes and alias parity.
- Add smoke coverage for default/verbose file-output modes.

Acceptance:
- CI covers at least 3.11, 3.12, 3.13 for tests and quality checks.
- All existing and new tests pass.

### PR 6: Security and Release Guardrails
- Add scheduled `pip-audit` CI workflow.
- Document current `biopython` CVE status and tracking note until fixed upstream.
- Add release checklist in docs.

Acceptance:
- Scheduled security job exists and reports status.
- Release docs cover lint/type/test/audit steps.

## Backlog (Post-Plan)
- Migrate from `flake8` to `ruff` fully for lint performance.
- Add benchmark regression guard for key commands.
- Introduce property-based tests for parser/formatter services.
