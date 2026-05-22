"""Tests for the IO layer: mutation parsing, structure inspection, output writing."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from stability_predictor_core.io import (
    MutationInputError,
    inspect_structure,
    load_mutations,
    sha256_file,
    write_predictions_json,
    write_result_json,
)


# ---------------------------------------------------------------------------
# load_mutations
# ---------------------------------------------------------------------------

def test_load_mutations_basic(tmp_path: Path) -> None:
    path = tmp_path / "muts.json"
    path.write_text(json.dumps({
        "chain": "A",
        "mutations": [
            {"position": 328, "wt": "A", "mt": "S"},
            {"position": 508, "wt": "F", "mt": "A"},
        ],
    }))
    muts = load_mutations(path)
    assert len(muts) == 2
    assert muts[0].chain == "A"
    assert muts[0].position == 328
    assert muts[1].wt == "F"
    assert muts[1].mt == "A"


def test_load_mutations_per_entry_chain_override(tmp_path: Path) -> None:
    path = tmp_path / "muts.json"
    path.write_text(json.dumps({
        "chain": "A",
        "mutations": [
            {"position": 1, "wt": "M", "mt": "L"},
            {"position": 5, "wt": "V", "mt": "I", "chain": "B"},
        ],
    }))
    muts = load_mutations(path)
    assert muts[0].chain == "A"
    assert muts[1].chain == "B"


def test_load_mutations_uppercases_aa(tmp_path: Path) -> None:
    path = tmp_path / "muts.json"
    path.write_text(json.dumps({
        "chain": "A",
        "mutations": [{"position": 1, "wt": "a", "mt": "v"}],
    }))
    muts = load_mutations(path)
    assert muts[0].wt == "A"
    assert muts[0].mt == "V"


@pytest.mark.parametrize("payload", [
    "{",  # malformed JSON
    json.dumps([]),  # not a dict at top
    json.dumps({"chain": "A"}),  # no mutations key
    json.dumps({"chain": "A", "mutations": []}),  # empty list
    json.dumps({"chain": "A", "mutations": [{"wt": "F", "mt": "A"}]}),  # missing position
    json.dumps({"mutations": [{"position": 1, "wt": "A", "mt": "V"}]}),  # no chain anywhere
])
def test_load_mutations_invalid_inputs_raise(tmp_path: Path, payload: str) -> None:
    path = tmp_path / "muts.json"
    path.write_text(payload)
    with pytest.raises(MutationInputError):
        load_mutations(path)


def test_load_mutations_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(MutationInputError):
        load_mutations(tmp_path / "does_not_exist.json")


# ---------------------------------------------------------------------------
# inspect_structure
# ---------------------------------------------------------------------------

def test_inspect_structure_minimal_pdb(tmp_path: Path) -> None:
    """The fallback parser handles CA-only stub PDBs (as used in tests)."""
    pdb = tmp_path / "stub.pdb"
    pdb.write_text(
        "HEADER    test\n"
        "ATOM      1  CA  ALA A 328       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  PHE B 508       3.800   0.000   0.000  1.00  0.00           C  \n"
        "END\n"
    )
    summary = inspect_structure(pdb)
    assert summary.residue_count == 2
    assert set(summary.chains) == {"A", "B"}
    text = summary.human_readable()
    assert "2 residues" in text
    assert "A" in text and "B" in text


def test_inspect_structure_human_readable_focus(tmp_path: Path) -> None:
    pdb = tmp_path / "stub.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA B   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "END\n"
    )
    summary = inspect_structure(pdb)
    text = summary.human_readable(focus_chain="A")
    assert "A" in text


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def test_sha256_file_stable(tmp_path: Path) -> None:
    f = tmp_path / "x.txt"
    f.write_text("hello")
    h1 = sha256_file(f)
    h2 = sha256_file(f)
    assert h1 == h2
    assert len(h1) == 64


def test_write_result_json_shape(tmp_path: Path) -> None:
    from stability_predictor_core.consensus import (
        ConsensusResult, Confidence, Direction,
    )
    from stability_predictor_core.methods.base import Mutation, StabilityPrediction

    mut = Mutation("A", 1, "A", "V")
    pred = StabilityPrediction(
        mutation=mut, ddg=0.1, confidence=0.7, method="rasp", method_version="x",
    )
    cr = ConsensusResult(
        mutation=mut, consensus_ddg=0.1, direction=Direction.NEUTRAL,
        confidence=Confidence.MEDIUM, n_methods_succeeded=1, n_methods_attempted=1,
        methods_agreed_direction=True, methods_agreed_magnitude=True,
        per_method=[pred], flags=[],
    )
    write_result_json(
        tmp_path, consensus_results=[cr],
        structure_path=Path("fake.pdb"),
        methods_requested=["rasp"], methods_available=["rasp"],
    )
    data = json.loads((tmp_path / "result.json").read_text())
    assert data["skill"] == "stability-predictor"
    assert data["input"]["n_mutations"] == 1
    assert data["summary"]["n_neutral"] == 1
    assert data["methods"]["skipped"] == []


def test_write_predictions_json_shape(tmp_path: Path) -> None:
    from stability_predictor_core.consensus import (
        ConsensusResult, Confidence, Direction,
    )
    from stability_predictor_core.methods.base import Mutation, StabilityPrediction

    mut = Mutation("A", 1, "A", "V")
    pred = StabilityPrediction(
        mutation=mut, ddg=0.1, confidence=0.7, method="rasp", method_version="x",
    )
    cr = ConsensusResult(
        mutation=mut, consensus_ddg=0.1, direction=Direction.NEUTRAL,
        confidence=Confidence.MEDIUM, n_methods_succeeded=1, n_methods_attempted=1,
        methods_agreed_direction=True, methods_agreed_magnitude=True,
        per_method=[pred], flags=[],
    )
    write_predictions_json(tmp_path, [cr])
    data = json.loads((tmp_path / "predictions.json").read_text())
    assert len(data["predictions"]) == 1
    assert data["predictions"][0]["mutation"] == "A:A1V"
