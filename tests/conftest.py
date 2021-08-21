from pathlib import Path

here = Path(__file__)


def pytest_configure(config):
    config.addinivalue_line("markers", "integration: mark as integration test")
