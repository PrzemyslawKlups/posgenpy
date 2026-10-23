import json
import os
from pathlib import Path

import pytest

CONFIG_ENV = "POSGENPY_PARITY_CONFIG"
DEFAULT_CONFIG = Path(__file__).parent / "parity_datasets.json"


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: needs the reference posgen datasets (see tests/README.md)"
    )


def load_parity_datasets():
    """Read the parity dataset config, or return [] if it isn't set up.

    Point POSGENPY_PARITY_CONFIG at a JSON file, or drop one at
    tests/parity_datasets.json (gitignored - it holds absolute paths to data
    that does not live in this repo). See parity_datasets.example.json.
    """
    path = Path(os.environ.get(CONFIG_ENV, DEFAULT_CONFIG))

    if not path.exists():
        return []

    with open(path) as f:
        datasets = json.load(f)

    return [d for d in datasets if Path(d["pos"]).exists()]


def pytest_generate_tests(metafunc):
    if "parity_dataset" in metafunc.fixturenames:
        datasets = load_parity_datasets()
        if not datasets:
            metafunc.parametrize("parity_dataset", [
                pytest.param(
                    None,
                    marks=pytest.mark.skip(
                        reason=f"no reference datasets; set ${CONFIG_ENV} "
                               f"or create {DEFAULT_CONFIG.name}"
                    ),
                )
            ])
        else:
            metafunc.parametrize(
                "parity_dataset", datasets, ids=[d["name"] for d in datasets]
            )
