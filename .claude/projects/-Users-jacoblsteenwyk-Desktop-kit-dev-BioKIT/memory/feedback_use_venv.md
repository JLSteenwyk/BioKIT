---
name: Use local venv for installs
description: Always install packages and run tools in the local venv virtual environment
type: feedback
---

Install packages in the local virtual environment named `venv`, not system Python or conda.

**Why:** User preference for keeping dependencies isolated in the project-local venv.

**How to apply:** Activate with `source venv/bin/activate` or use `venv/bin/python` / `venv/bin/pip` directly before running tests, linters, or installing packages.
