---
name: BioKIT Overview
description: Core architecture and structure of the BioKIT bioinformatics toolkit
type: project
---

BioKIT is a Python CLI toolkit for molecular sequence data analysis, published in Genetics (2022). v1.1.5, Python 3.11-3.13, MIT license.

**Architecture:** Service-based — `biokit.py` dispatches CLI commands to service classes under `biokit/services/` (alignment, coding_sequences, fastq, genome, text). All services inherit from `BaseService` which provides structured output (TSV/JSON/YAML).

**Key paths:**
- CLI dispatcher: `biokit/biokit.py`
- Services: `biokit/services/{category}/{command}.py`
- Helpers: `biokit/helpers/` (file I/O, character frequencies)
- Tests: `tests/unit/` and `tests/integration/` with 145 sample fixtures
- Dev scripts: `scripts/` (benchmarking, private API guard)

**Dev workflow:** `make dev.setup`, `make quality` (lint + typecheck + guard + test.fast), `make test` for full suite. CI via GitHub Actions with quality, fast-test (3.11-3.13 matrix), and full-coverage jobs.

**Why:** Understanding the architecture helps me navigate the codebase quickly and make changes consistent with existing patterns.

**How to apply:** Follow the service pattern (BaseService subclass with `process_args()` and `run()`) when adding commands. Use `make quality` before suggesting commits.
