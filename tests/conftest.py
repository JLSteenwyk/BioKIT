import sys
from pathlib import Path

import pytest

here = Path(__file__)


def pytest_configure(config):
    config.addinivalue_line("markers", "integration: mark as integration test")


@pytest.fixture(autouse=True)
def reset_sys_argv(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["biokit"])
