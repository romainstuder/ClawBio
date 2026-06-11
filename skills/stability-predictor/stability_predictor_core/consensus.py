"""Consensus logic: aggregate predictions across multiple methods.

A single mutation may be predicted by RaSP, ThermoMPNN, and FoldX. This module:
  1. Reduces them to a consensus ΔΔG (mean of successful predictions).
  2. Classifies direction: destabilizing / neutral / stabilizing.
  3. Assigns a confidence tier based on inter-method agreement.
  4. Flags disagreements that warrant manual review.

Thresholds are conservative defaults; users can override via CLI flags.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from statistics import mean

from .methods.base import Mutation, StabilityPrediction

# Default thresholds (kcal/mol). Tunable but rarely should be.
NEUTRAL_BAND_KCAL_MOL = 1.0
"""ΔΔG within ±1.0 kcal/mol is considered neutral.

Roughly one hydrogen bond's worth of energy. Standard threshold in the literature
(e.g., Tokuriki & Tawfik 2009 Curr Opin Struct Biol).
"""

AGREEMENT_TOLERANCE_KCAL_MOL = 1.0
"""Methods are said to 'agree on magnitude' if they're within this range.

Loose threshold; stricter agreement (0.5) would flag most multi-method runs.
"""


class Direction(str, Enum):
    DESTABILIZING = "destabilizing"
    NEUTRAL = "neutral"
    STABILIZING = "stabilizing"
    UNKNOWN = "unknown"


class Confidence(str, Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    NONE = "none"


@dataclass
class ConsensusResult:
    """Aggregated prediction for one mutation across multiple methods."""

    mutation: Mutation
    consensus_ddg: float | None
    direction: Direction
    confidence: Confidence
    n_methods_succeeded: int
    n_methods_attempted: int
    methods_agreed_direction: bool
    methods_agreed_magnitude: bool
    per_method: list[StabilityPrediction]
    flags: list[str]

    def to_dict(self) -> dict:
        return {
            "mutation": str(self.mutation),
            "consensus_ddg": self.consensus_ddg,
            "direction": self.direction.value,
            "confidence": self.confidence.value,
            "n_methods_succeeded": self.n_methods_succeeded,
            "n_methods_attempted": self.n_methods_attempted,
            "methods_agreed_direction": self.methods_agreed_direction,
            "methods_agreed_magnitude": self.methods_agreed_magnitude,
            "per_method": [p.to_dict() for p in self.per_method],
            "flags": self.flags,
        }


def classify_direction(ddg: float, neutral_band: float = NEUTRAL_BAND_KCAL_MOL) -> Direction:
    """Classify a ΔΔG value into destabilizing / neutral / stabilizing."""
    if ddg > neutral_band:
        return Direction.DESTABILIZING
    if ddg < -neutral_band:
        return Direction.STABILIZING
    return Direction.NEUTRAL


def aggregate(
    predictions_per_mutation: list[list[StabilityPrediction]],
    neutral_band: float = NEUTRAL_BAND_KCAL_MOL,
    agreement_tolerance: float = AGREEMENT_TOLERANCE_KCAL_MOL,
) -> list[ConsensusResult]:
    """Aggregate per-method predictions into per-mutation consensus.

    Args:
        predictions_per_mutation: One inner list per mutation, containing
            predictions from each requested method.
        neutral_band: ΔΔG magnitude below which a mutation is called neutral.
        agreement_tolerance: kcal/mol within which methods are considered
            to agree on magnitude.

    Returns:
        One ConsensusResult per input mutation, in input order.
    """
    return [
        _aggregate_one(preds, neutral_band, agreement_tolerance)
        for preds in predictions_per_mutation
    ]


def _aggregate_one(
    predictions: list[StabilityPrediction],
    neutral_band: float,
    agreement_tolerance: float,
) -> ConsensusResult:
    if not predictions:
        raise ValueError("Cannot aggregate empty prediction list")

    mutation = predictions[0].mutation
    if not all(p.mutation == mutation for p in predictions):
        raise ValueError("All predictions in a group must be for the same mutation")

    successful = [p for p in predictions if p.succeeded]
    n_succeeded = len(successful)
    n_attempted = len(predictions)
    flags: list[str] = []

    if n_succeeded == 0:
        flags.extend(_collect_error_flags(predictions))
        return ConsensusResult(
            mutation=mutation,
            consensus_ddg=None,
            direction=Direction.UNKNOWN,
            confidence=Confidence.NONE,
            n_methods_succeeded=0,
            n_methods_attempted=n_attempted,
            methods_agreed_direction=False,
            methods_agreed_magnitude=False,
            per_method=predictions,
            flags=flags,
        )

    ddg_values: list[float] = [p.ddg for p in successful if p.ddg is not None]
    consensus = mean(ddg_values)
    direction = classify_direction(consensus, neutral_band)

    directions = {classify_direction(d, neutral_band) for d in ddg_values}
    agreed_direction = len(directions) == 1
    agreed_magnitude = _agreed_on_magnitude(ddg_values, agreement_tolerance)

    confidence = _compute_confidence(
        n_succeeded=n_succeeded,
        n_attempted=n_attempted,
        agreed_direction=agreed_direction,
        agreed_magnitude=agreed_magnitude,
    )

    if not agreed_direction and n_succeeded >= 2:
        flags.append("Methods disagree on direction; review individual predictions")
    if not agreed_magnitude and n_succeeded >= 2:
        spread = max(ddg_values) - min(ddg_values)
        flags.append(f"Method spread {spread:.2f} kcal/mol exceeds tolerance")
    if n_succeeded < n_attempted:
        flags.append(f"{n_attempted - n_succeeded} of {n_attempted} methods failed")

    return ConsensusResult(
        mutation=mutation,
        consensus_ddg=consensus,
        direction=direction,
        confidence=confidence,
        n_methods_succeeded=n_succeeded,
        n_methods_attempted=n_attempted,
        methods_agreed_direction=agreed_direction,
        methods_agreed_magnitude=agreed_magnitude,
        per_method=predictions,
        flags=flags,
    )


def _agreed_on_magnitude(values: list[float], tolerance: float) -> bool:
    """All values within `tolerance` of each other (max - min ≤ tolerance)."""
    if len(values) < 2:
        return True
    return (max(values) - min(values)) <= tolerance


def _compute_confidence(
    *,
    n_succeeded: int,
    n_attempted: int,
    agreed_direction: bool,
    agreed_magnitude: bool,
) -> Confidence:
    """Confidence rubric.

    HIGH:   ≥2 methods succeed AND agree on both direction and magnitude
    MEDIUM: ≥2 methods succeed AND agree on direction (magnitude may vary), OR
            single method succeeded but no failures
    LOW:    Single method succeeded with at least one failure, OR
            multiple methods succeeded but disagree on direction
    NONE:   All methods failed
    """
    if n_succeeded == 0:
        return Confidence.NONE
    if n_succeeded >= 2 and agreed_direction and agreed_magnitude:
        return Confidence.HIGH
    if n_succeeded >= 2 and agreed_direction:
        return Confidence.MEDIUM
    if n_succeeded == 1 and n_attempted == 1:
        # Only one method requested and it worked — medium, not high.
        # User opted not to cross-check; we don't penalise that but we don't reward it either.
        return Confidence.MEDIUM
    if n_succeeded == 1:
        return Confidence.LOW
    # ≥2 succeeded but disagree on direction
    return Confidence.LOW


def _collect_error_flags(predictions: list[StabilityPrediction]) -> list[str]:
    """Build human-readable flags for a fully-failed mutation."""
    flags: list[str] = ["All methods failed for this mutation"]
    for p in predictions:
        if p.error:
            flags.append(f"  {p.method}: {p.error}")
    return flags
