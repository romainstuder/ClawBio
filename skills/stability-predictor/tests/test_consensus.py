"""Tests for the consensus aggregation logic.

Mirrors the original sanity_test.py but as pytest test functions.
"""

from __future__ import annotations

import pytest

from stability_predictor_core.consensus import (
    Confidence,
    Direction,
    aggregate,
    classify_direction,
)
from stability_predictor_core.methods.base import Mutation, StabilityPrediction


def _pred(method: str, mutation: Mutation, ddg: float | None,
          error: str | None = None) -> StabilityPrediction:
    return StabilityPrediction(
        mutation=mutation,
        ddg=ddg,
        confidence=0.7 if ddg is not None else None,
        method=method,
        method_version="test",
        error=error,
    )


def test_classify_direction_thresholds():
    assert classify_direction(2.0) == Direction.DESTABILIZING
    assert classify_direction(0.5) == Direction.NEUTRAL
    assert classify_direction(-2.0) == Direction.STABILIZING
    # exactly on the boundary is neutral (band is inclusive)
    assert classify_direction(1.0) == Direction.NEUTRAL


def test_three_methods_agree_destabilizing():
    """F508del-like case: all three methods agree on destabilising."""
    mut = Mutation(chain="A", position=508, wt="F", mt="A")
    preds = [
        _pred("rasp", mut, 3.2),
        _pred("thermompnn", mut, 2.9),
        _pred("foldx", mut, 3.4),
    ]
    [result] = aggregate([preds])
    assert result.direction == Direction.DESTABILIZING
    assert result.confidence == Confidence.HIGH
    assert result.methods_agreed_direction
    assert result.methods_agreed_magnitude
    assert result.flags == []
    assert abs(result.consensus_ddg - 3.17) < 0.01


def test_direction_disagreement_lowers_confidence():
    """RaSP says stabilizing, FoldX says destabilizing → low confidence + flag."""
    mut = Mutation(chain="A", position=100, wt="V", mt="L")
    preds = [_pred("rasp", mut, -1.5), _pred("foldx", mut, 1.8)]
    [result] = aggregate([preds])
    assert not result.methods_agreed_direction
    assert result.confidence == Confidence.LOW
    assert any("disagree on direction" in f for f in result.flags)


def test_magnitude_disagreement_still_medium():
    """Same direction, spread > tolerance → MEDIUM with a spread flag."""
    mut = Mutation(chain="A", position=200, wt="G", mt="W")
    preds = [_pred("rasp", mut, 1.2), _pred("foldx", mut, 3.5)]
    [result] = aggregate([preds])
    assert result.methods_agreed_direction
    assert not result.methods_agreed_magnitude
    assert result.confidence == Confidence.MEDIUM
    assert any("spread" in f for f in result.flags)


def test_one_method_fails_one_succeeds():
    mut = Mutation(chain="A", position=300, wt="K", mt="E")
    preds = [
        _pred("rasp", mut, None, error="model weights not found"),
        _pred("foldx", mut, 0.5),
    ]
    [result] = aggregate([preds])
    assert result.n_methods_succeeded == 1
    assert result.n_methods_attempted == 2
    assert result.confidence == Confidence.LOW
    assert result.direction == Direction.NEUTRAL


def test_all_methods_fail():
    mut = Mutation(chain="A", position=400, wt="P", mt="A")
    preds = [
        _pred("rasp", mut, None, error="timeout"),
        _pred("foldx", mut, None, error="binary not found"),
    ]
    [result] = aggregate([preds])
    assert result.consensus_ddg is None
    assert result.direction == Direction.UNKNOWN
    assert result.confidence == Confidence.NONE
    assert any("All methods failed" in f for f in result.flags)


def test_single_method_run_is_medium_not_high():
    """User chose only one method; we don't promote to HIGH without cross-check."""
    mut = Mutation(chain="A", position=50, wt="A", mt="V")
    [result] = aggregate([[_pred("rasp", mut, 0.3)]])
    assert result.confidence == Confidence.MEDIUM
    assert result.direction == Direction.NEUTRAL


def test_mismatched_mutation_raises():
    """All predictions in a group must be for the same mutation."""
    m1 = Mutation("A", 1, "A", "V")
    m2 = Mutation("A", 2, "L", "I")
    with pytest.raises(ValueError):
        aggregate([[_pred("rasp", m1, 0.1), _pred("foldx", m2, 0.2)]])
