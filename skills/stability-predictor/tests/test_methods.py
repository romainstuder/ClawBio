"""Tests for the method wrappers — availability, install instructions, and
graceful failure when the real backend isn't installed.

Real-model tests live in `test_validation.py` and are marked `@pytest.mark.slow`.
"""

from __future__ import annotations

import json
from pathlib import Path

from stability_predictor_core.methods.base import (
    Mutation,
    StabilityMethod,
)
from stability_predictor_core.methods.foldx import FoldX
from stability_predictor_core.methods.precomputed import Precomputed
from stability_predictor_core.methods.rasp import RaSP
from stability_predictor_core.methods.thermompnn import ThermoMPNN


# ---------------------------------------------------------------------------
# Interface conformance
# ---------------------------------------------------------------------------

def test_all_methods_implement_interface():
    """Each wrapper exposes name, version, is_available, install_instructions."""
    for method_cls in (RaSP, ThermoMPNN, FoldX):
        m: StabilityMethod = method_cls()
        assert isinstance(m.name, str) and m.name
        assert isinstance(m.version, str) and m.version
        assert isinstance(m.is_available(), bool)
        instructions = m.install_instructions()
        assert isinstance(instructions, str) and instructions.strip()


def test_install_instructions_mention_install_step():
    """Each method's install instructions tell the user how to install it."""
    assert "pip install" in RaSP().install_instructions().lower()
    assert "pip install" in ThermoMPNN().install_instructions().lower()
    assert "foldx" in FoldX().install_instructions().lower()


# ---------------------------------------------------------------------------
# Graceful unavailability
# ---------------------------------------------------------------------------

def test_rasp_without_backend_returns_failed_predictions(tmp_path: Path) -> None:
    """If the rasp-predictor package isn't installed, calls return error predictions."""
    method = RaSP()
    if method.is_available():
        return  # In a real-env CI we'd skip; here we just no-op.
    pdb = tmp_path / "stub.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "END\n"
    )
    preds = method.predict(pdb, [Mutation("A", 1, "A", "V")])
    assert len(preds) == 1
    assert preds[0].error is not None
    assert not preds[0].succeeded


def test_thermompnn_without_backend_returns_failed_predictions(tmp_path: Path) -> None:
    method = ThermoMPNN()
    if method.is_available():
        return
    pdb = tmp_path / "stub.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "END\n"
    )
    preds = method.predict(pdb, [Mutation("A", 1, "A", "V")])
    assert preds[0].error is not None


def test_foldx_without_binary_returns_failed_predictions(tmp_path: Path) -> None:
    method = FoldX()
    if method.is_available():
        return
    pdb = tmp_path / "stub.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "END\n"
    )
    preds = method.predict(pdb, [Mutation("A", 1, "A", "V")])
    assert preds[0].error is not None
    assert "foldx" in (preds[0].error or "").lower()


# ---------------------------------------------------------------------------
# Precomputed (used for demo)
# ---------------------------------------------------------------------------

def test_precomputed_loads_and_predicts(tmp_path: Path) -> None:
    cache = tmp_path / "cache.json"
    cache.write_text(json.dumps({
        "method": "rasp",
        "method_version": "1.2.0 (cached)",
        "predictions": [
            {"chain": "A", "position": 1, "wt": "A", "mt": "V",
             "ddg": 0.5, "notes": "near-neutral"},
        ],
    }))
    m = Precomputed(cache)
    assert m.is_available()
    assert m.name == "rasp"
    preds = m.predict(Path("ignored.pdb"), [Mutation("A", 1, "A", "V")])
    assert preds[0].succeeded
    assert preds[0].ddg == 0.5


def test_precomputed_missing_cache_is_unavailable(tmp_path: Path) -> None:
    m = Precomputed(tmp_path / "nope.json")
    assert not m.is_available()
    assert "not found" in m.install_instructions().lower()


def test_precomputed_unknown_mutation_returns_error(tmp_path: Path) -> None:
    cache = tmp_path / "cache.json"
    cache.write_text(json.dumps({
        "method": "rasp", "method_version": "x",
        "predictions": [
            {"chain": "A", "position": 1, "wt": "A", "mt": "V", "ddg": 0.5},
        ],
    }))
    m = Precomputed(cache)
    preds = m.predict(Path("ignored.pdb"), [Mutation("A", 99, "L", "I")])
    assert preds[0].error is not None
    assert not preds[0].succeeded


def test_precomputed_label_overrides_name(tmp_path: Path) -> None:
    cache = tmp_path / "cache.json"
    cache.write_text(json.dumps({
        "method": "rasp", "method_version": "x",
        "predictions": [{"chain": "A", "position": 1, "wt": "A", "mt": "V", "ddg": 0.0}],
    }))
    m = Precomputed(cache, label="foldx")
    assert m.name == "foldx"
