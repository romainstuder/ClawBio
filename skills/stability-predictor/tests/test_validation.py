"""Tier-3 benchmark tests — slow, marked, and skipped by default.

These exercise the actual prediction backends (RaSP, ThermoMPNN, FoldX) and
measure Spearman correlation against an experimental ΔΔG benchmark
(ProTherm/FireProtDB subset). They require:

  * The relevant package or binary installed
  * Network access and ~50-200 MB of model weights cached
  * Minutes of runtime

Run explicitly with:

    pytest skills/stability-predictor/tests/test_validation.py -m slow

Expected thresholds, taken from the upstream publications, are:
  RaSP        Spearman >= 0.55
  ThermoMPNN  Spearman >= 0.60
  FoldX       Spearman >= 0.50

Implementation note: the benchmark dataset is not vendored in this repo to
keep the checkout light. When run for real the test fetches a small curated
subset and caches it under ~/.cache/clawbio/protherm_subset.tsv.
"""

from __future__ import annotations

import pytest

from stability_predictor_core.methods.foldx import FoldX
from stability_predictor_core.methods.rasp import RaSP
from stability_predictor_core.methods.thermompnn import ThermoMPNN

pytestmark = pytest.mark.slow


@pytest.mark.skipif(not RaSP().is_available(), reason="RaSP not installed")
def test_rasp_protherm_correlation() -> None:
    pytest.skip("ProTherm benchmark fetch not implemented in v1 scaffold")


@pytest.mark.skipif(not ThermoMPNN().is_available(), reason="ThermoMPNN not installed")
def test_thermompnn_protherm_correlation() -> None:
    pytest.skip("ProTherm benchmark fetch not implemented in v1 scaffold")


@pytest.mark.skipif(not FoldX().is_available(), reason="FoldX binary not on PATH")
def test_foldx_protherm_correlation() -> None:
    pytest.skip("ProTherm benchmark fetch not implemented in v1 scaffold")
