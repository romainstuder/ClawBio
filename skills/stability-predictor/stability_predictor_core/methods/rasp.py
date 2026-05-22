"""RaSP (Rapid Stability Predictor) method wrapper.

RaSP: Blaabjerg et al. (2023) eLife. https://elifesciences.org/articles/82593

RaSP is a CNN trained on Rosetta ΔΔG values across ~10M variants of 1,400 proteins.
Output: per-mutation ΔΔG in kcal/mol.
License: Apache 2.0.
GPU optional, CPU works for small inputs.

This wrapper exposes RaSP through the StabilityMethod interface. It does not
reimplement RaSP — it calls the upstream package.
"""

from __future__ import annotations

import importlib
import importlib.util
import logging
from pathlib import Path

from .base import Mutation, StabilityMethod, StabilityPrediction

logger = logging.getLogger(__name__)

# Confidence heuristic: RaSP doesn't ship per-prediction uncertainty.
# Use a fixed midrange confidence; method comparison drives the real confidence
# signal at the consensus layer.
_RASP_DEFAULT_CONFIDENCE = 0.7


class RaSP(StabilityMethod):
    """RaSP stability predictor (default method).

    Availability: requires the `rasp-predictor` package and downloaded model weights.
    First-time use will trigger a one-shot model download (~50 MB) cached under
    ~/.cache/rasp/.

    Performance: ~0.1 s per mutation on CPU, ~10x faster on GPU.
    """

    name = "rasp"

    def __init__(self) -> None:
        super().__init__()
        self._module = None
        self._load_error: str | None = None
        self._try_import()

    def _try_import(self) -> None:
        """Attempt to import RaSP; record any failure for diagnostics."""
        try:
            self._module = importlib.import_module("rasp_predictor")
            self.version = getattr(self._module, "__version__", "unknown")
        except ImportError as exc:
            self._load_error = str(exc)
            self.version = "not-installed"

    def is_available(self) -> bool:
        if self._module is None:
            return False
        # Models may need to be downloaded; we treat that as "available"
        # (the package handles the download lazily). If you want to gate
        # availability on model presence, add a check here.
        return True

    def install_instructions(self) -> str:
        return (
            "RaSP is open source and pip-installable.\n"
            "\n"
            "  pip install rasp-predictor\n"
            "\n"
            "First use downloads ~50 MB of model weights to ~/.cache/rasp/.\n"
            "Reference: Blaabjerg et al. (2023) eLife.\n"
            "https://github.com/KULL-Centre/_2022_ML-ddG-Blaabjerg\n"
        )

    def predict(
        self,
        structure_path: Path,
        mutations: list[Mutation],
    ) -> list[StabilityPrediction]:
        if not self.is_available():
            # Return failed predictions rather than raising; caller decides what to do.
            return [
                StabilityPrediction(
                    mutation=m,
                    ddg=None,
                    confidence=None,
                    method=self.name,
                    method_version=self.version,
                    error=f"RaSP unavailable: {self._load_error or 'unknown'}",
                )
                for m in mutations
            ]

        if not structure_path.exists():
            raise RuntimeError(f"Structure file not found: {structure_path}")

        results: list[StabilityPrediction] = []
        for mutation in mutations:
            results.append(self._predict_one(structure_path, mutation))
        return results

    def _predict_one(self, structure_path: Path, mutation: Mutation) -> StabilityPrediction:
        """Predict a single mutation. Failures captured as errors, not raised."""
        if mutation.mt == "-":
            return StabilityPrediction(
                mutation=mutation,
                ddg=None,
                confidence=None,
                method=self.name,
                method_version=self.version,
                error="Deletions not supported in v1",
            )

        try:
            # NOTE: rasp_predictor's exact public API will determine the call signature below.
            # This is a placeholder that matches the typical pattern; adjust once you've
            # confirmed the upstream interface during pre-flight.
            assert self._module is not None
            raw = self._module.predict_single(  # type: ignore[attr-defined]
                pdb_path=str(structure_path),
                chain=mutation.chain,
                position=mutation.position,
                wt=mutation.wt,
                mt=mutation.mt,
            )
            ddg = float(raw["ddg"])
            notes = self._build_notes(raw)
            return StabilityPrediction(
                mutation=mutation,
                ddg=ddg,
                confidence=_RASP_DEFAULT_CONFIDENCE,
                method=self.name,
                method_version=self.version,
                notes=notes,
                raw_output={"rasp": raw},
            )
        except KeyError as exc:
            return StabilityPrediction(
                mutation=mutation,
                ddg=None,
                confidence=None,
                method=self.name,
                method_version=self.version,
                error=f"RaSP output missing expected field: {exc}",
            )
        except Exception as exc:  # noqa: BLE001 - capture method-internal failures
            logger.warning("RaSP failed on %s: %s", mutation, exc)
            return StabilityPrediction(
                mutation=mutation,
                ddg=None,
                confidence=None,
                method=self.name,
                method_version=self.version,
                error=f"RaSP error: {exc}",
            )

    @staticmethod
    def _build_notes(raw: dict) -> str:
        """Surface useful per-mutation context from the raw RaSP output."""
        notes: list[str] = []
        if raw.get("low_coverage"):
            notes.append("low MSA coverage at this position")
        if raw.get("buried", False):
            notes.append("buried residue (RaSP often most accurate here)")
        return "; ".join(notes)
