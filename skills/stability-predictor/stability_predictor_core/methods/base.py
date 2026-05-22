"""Abstract base class for stability prediction methods.

All methods (RaSP, ThermoMPNN, FoldX, future additions) implement this interface.
The CLI dispatches to one or more methods via this contract.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class Mutation:
    """A single point mutation specification.

    Attributes:
        chain: Chain identifier (e.g., "A"). Required.
        position: 1-based residue number in the chain.
        wt: Wild-type one-letter amino acid code.
        mt: Mutant one-letter amino acid code. Use "-" for deletion (not supported in v1).
    """

    chain: str
    position: int
    wt: str
    mt: str

    def __str__(self) -> str:
        return f"{self.chain}:{self.wt}{self.position}{self.mt}"


@dataclass
class StabilityPrediction:
    """Result of a single ΔΔG prediction.

    Attributes:
        mutation: The mutation predicted.
        ddg: Folding free energy change in kcal/mol.
                Positive = destabilizing, negative = stabilizing.
                None if prediction failed.
        confidence: Method-specific confidence in [0, 1]. None if not applicable.
        method: Method name (e.g., "rasp").
        method_version: Method version string for reproducibility.
        notes: Human-readable notes (e.g., "low coverage at this site").
        error: Error message if prediction failed; None on success.
        raw_output: Method-specific raw data, for reproducibility/debugging.
    """

    mutation: Mutation
    ddg: float | None
    confidence: float | None
    method: str
    method_version: str
    notes: str = ""
    error: str | None = None
    raw_output: dict[str, Any] = field(default_factory=dict)

    @property
    def succeeded(self) -> bool:
        return self.ddg is not None and self.error is None

    def to_dict(self) -> dict[str, Any]:
        """Serialize for JSON output. Excludes raw_output (often large)."""
        return {
            "mutation": str(self.mutation),
            "chain": self.mutation.chain,
            "position": self.mutation.position,
            "wt": self.mutation.wt,
            "mt": self.mutation.mt,
            "ddg": self.ddg,
            "confidence": self.confidence,
            "method": self.method,
            "method_version": self.method_version,
            "notes": self.notes,
            "error": self.error,
        }


class StabilityMethod(ABC):
    """Base class for all stability prediction methods.

    Concrete subclasses must implement:
      - name (class attr): short identifier used in CLI (e.g., "rasp")
      - version (instance attr): version string after init
      - is_available(): can this method run on this machine?
      - install_instructions(): user-facing setup help
      - predict(): the actual prediction logic
    """

    name: str = "abstract"

    def __init__(self) -> None:
        self.version: str = "unknown"

    @abstractmethod
    def is_available(self) -> bool:
        """Return True if all dependencies (binaries, models, weights) are present.

        Must not raise. Used by the CLI to decide whether to attempt prediction
        and what install instructions to show.
        """

    @abstractmethod
    def install_instructions(self) -> str:
        """Return human-readable setup instructions.

        Shown to the user when they request this method but it's unavailable.
        Should include exact commands, URLs, and a 'verify' step.
        """

    @abstractmethod
    def predict(
        self,
        structure_path: Path,
        mutations: list[Mutation],
    ) -> list[StabilityPrediction]:
        """Predict ΔΔG for each mutation on the given structure.

        Args:
            structure_path: Path to PDB or CIF file.
            mutations: List of mutations to predict. Order preserved in output.

        Returns:
            One StabilityPrediction per input mutation, in the same order.
            Per-mutation failures populate `error` but do not raise.

        Raises:
            RuntimeError: Only for catastrophic failures (e.g., structure unreadable).
                Per-mutation prediction errors must be captured in the result, not raised.
        """

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(version={self.version})"
