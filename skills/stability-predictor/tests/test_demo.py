"""End-to-end CLI smoke tests using the bundled demo data.

These do NOT require RaSP, ThermoMPNN, or FoldX to be installed — they run
the Precomputed cache method that ships with the demo data.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
CLI = SKILL_DIR / "stability_predictor.py"


def _run_demo(out_dir: Path, *extra: str) -> subprocess.CompletedProcess:
    """Invoke the CLI in a subprocess with --demo."""
    cmd = [sys.executable, str(CLI), "--demo", "--output", str(out_dir), *extra]
    return subprocess.run(cmd, capture_output=True, text=True)


@pytest.mark.parametrize("demo_set,n_mutations", [
    ("t4lysozyme", 3),
    ("p53", 1),
])
def test_demo_runs_and_writes_expected_files(
    tmp_path: Path, demo_set: str, n_mutations: int
) -> None:
    out = tmp_path / "out"
    result = _run_demo(out, "--demo-set", demo_set, "--method", "all")
    assert result.returncode == 0, f"CLI failed: {result.stderr}"

    expected = [
        "report.md",
        "result.json",
        "predictions.json",
        "reproducibility/commands.sh",
        "reproducibility/environment.yml",
        "reproducibility/input_checksum.txt",
        "reproducibility/output_checksum.txt",
    ]
    for relpath in expected:
        assert (out / relpath).exists(), f"missing output: {relpath}"

    result_json = json.loads((out / "result.json").read_text())
    assert result_json["skill"] == "stability-predictor"
    assert result_json["input"]["n_mutations"] == n_mutations
    # All three method names should appear when --method all + caches present
    assert set(result_json["methods"]["available"]) == {"rasp", "thermompnn", "foldx"}


def test_demo_report_has_required_sections(tmp_path: Path) -> None:
    out = tmp_path / "out"
    result = _run_demo(out, "--method", "all")
    assert result.returncode == 0

    report = (out / "report.md").read_text()
    for section in (
        "# Stability Prediction Report",
        "## Input",
        "## Summary",
        "## Interpretation Guide",
        "## Methods",
        "## References",
        "## Reproducibility",
        "## Disclaimer",
    ):
        assert section in report, f"missing section: {section}"


def test_demo_predictions_contain_all_mutations(tmp_path: Path) -> None:
    out = tmp_path / "out"
    result = _run_demo(out, "--demo-set", "t4lysozyme", "--method", "all")
    assert result.returncode == 0

    preds = json.loads((out / "predictions.json").read_text())
    assert len(preds["predictions"]) == 3
    mutations = {p["mutation"] for p in preds["predictions"]}
    # T4 lysozyme Matthews-lab benchmark set
    assert mutations == {"A:L99A", "A:T157I", "A:T26S"}


def test_demo_l99a_classified_destabilizing(tmp_path: Path) -> None:
    """L99A is the textbook destabilising mutation (~+5 kcal/mol)."""
    out = tmp_path / "out"
    result = _run_demo(out, "--demo-set", "t4lysozyme", "--method", "all")
    assert result.returncode == 0

    preds = json.loads((out / "predictions.json").read_text())
    by_mut = {p["mutation"]: p for p in preds["predictions"]}
    l99a = by_mut["A:L99A"]
    assert l99a["direction"] == "destabilizing"
    assert l99a["confidence"] == "high"
    assert l99a["consensus_ddg"] > 3.0


def test_demo_y220c_classified_destabilizing(tmp_path: Path) -> None:
    """Y220C is the canonical destabilising p53 cancer mutation."""
    out = tmp_path / "out"
    result = _run_demo(out, "--demo-set", "p53", "--method", "all")
    assert result.returncode == 0

    preds = json.loads((out / "predictions.json").read_text())
    [y220c] = preds["predictions"]
    assert y220c["mutation"] == "A:Y220C"
    assert y220c["direction"] == "destabilizing"
    assert y220c["consensus_ddg"] > 2.0


def test_demo_replay_script_is_portable(tmp_path: Path) -> None:
    """commands.sh should be built from parsed args, not sys.argv."""
    out = tmp_path / "out"
    result = _run_demo(out, "--method", "all")
    assert result.returncode == 0

    commands = (out / "reproducibility" / "commands.sh").read_text()
    assert commands.startswith("#!/usr/bin/env bash")
    assert "set -euo pipefail" in commands
    assert "--structure" in commands
    assert "--mutations" in commands
    assert "--method all" in commands


def test_demo_checksums_cover_all_outputs(tmp_path: Path) -> None:
    out = tmp_path / "out"
    result = _run_demo(out, "--method", "all")
    assert result.returncode == 0

    checksums = (out / "reproducibility" / "output_checksum.txt").read_text()
    for f in ("report.md", "result.json", "predictions.json"):
        assert f in checksums


def test_user_mode_fails_cleanly_without_inputs(tmp_path: Path) -> None:
    """Without --demo, --structure and --mutations are mandatory."""
    cmd = [sys.executable, str(CLI), "--output", str(tmp_path / "out")]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0
    assert "structure" in result.stderr.lower() or "demo" in result.stderr.lower()
