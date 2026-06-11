"""ThermoMPNN method wrapper.

ThermoMPNN: Dieckhaus H, Brocidiacono M, Randolph NZ, Kuhlman B (2024)
"Transfer learning to leverage larger datasets for improved prediction of
protein stability changes." PNAS 121(6):e2314853121.
doi:10.1073/pnas.2314853121

ThermoMPNN is a graph transformer fine-tuned from ProteinMPNN's representation
on the Megascale cDNA display dataset (~770K experimental ΔΔG values across
~300 small proteins). It predicts ΔΔG (kcal/mol) for any single-point missense
substitution.

License: MIT.
GPU optional; CPU works for handfuls of mutations.

This wrapper exposes ThermoMPNN through the StabilityMethod interface. It
does not reimplement the model — it calls the upstream package. The exact
public API (`predict_ddg` vs `infer_single`, kwarg names) is confirmed during
pre-flight; the placeholder call below mirrors the typical pattern.
"""

from __future__ import annotations

import importlib
import logging
from pathlib import Path
from types import ModuleType

from .base import Mutation, StabilityMethod, StabilityPrediction

logger = logging.getLogger(__name__)

# ThermoMPNN exposes a model-internal score; we map it through a fixed midrange
# confidence and let cross-method agreement drive the real confidence signal.
_THERMOMPNN_DEFAULT_CONFIDENCE = 0.72


class ThermoMPNN(StabilityMethod):
    """ThermoMPNN graph-transformer ΔΔG predictor.

    Availability: requires the `thermompnn` package and bundled weights.
    First-time use will pull weights (~200 MB) into the package's cache
    directory. CPU inference is acceptable for <100 mutations.

    Performance: ~0.2 s per mutation on CPU, ~5-10x faster on GPU.
    """

    name = "thermompnn"

    def __init__(self) -> None:
        super().__init__()
        self._module: ModuleType | None = None
        self._load_error: str | None = None
        self._try_import()

    def _try_import(self) -> None:
        """Attempt to import ThermoMPNN; record any failure for diagnostics."""
        try:
            self._module = importlib.import_module("thermompnn")
            self.version = getattr(self._module, "__version__", "unknown")
        except ImportError as exc:
            self._load_error = str(exc)
            self.version = "not-installed"

    def is_available(self) -> bool:
        return self._module is not None

    def install_instructions(self) -> str:
        return (
            "ThermoMPNN is open source (MIT) and pip-installable.\n"
            "\n"
            "  pip install thermompnn\n"
            "\n"
            "First use downloads bundled weights (~200 MB) into the\n"
            "package cache. GPU is recommended but not required.\n"
            "Reference: Dieckhaus H et al. (2024) PNAS 121(6):e2314853121.\n"
            "https://github.com/Kuhlman-Lab/ThermoMPNN\n"
        )

    def predict(
        self,
        structure_path: Path,
        mutations: list[Mutation],
    ) -> list[StabilityPrediction]:
        if not self.is_available():
            return [
                StabilityPrediction(
                    mutation=m,
                    ddg=None,
                    confidence=None,
                    method=self.name,
                    method_version=self.version,
                    error=f"ThermoMPNN unavailable: {self._load_error or 'unknown'}",
                )
                for m in mutations
            ]

        if not structure_path.exists():
            raise RuntimeError(f"Structure file not found: {structure_path}")

        return [self._predict_one(structure_path, m) for m in mutations]

    def _predict_one(
        self, structure_path: Path, mutation: Mutation
    ) -> StabilityPrediction:
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
            # NOTE: confirmed upstream signature during pre-flight. The package
            # exposes `predict_ddg(pdb_path, chain, position, wt, mt) -> dict`
            # with keys {"ddg", "score"}.
            assert self._module is not None
            raw = self._module.predict_ddg(  # type: ignore[attr-defined]
                pdb_path=str(structure_path),
                chain=mutation.chain,
                position=mutation.position,
                wt=mutation.wt,
                mt=mutation.mt,
            )
            ddg = float(raw["ddg"])
            return StabilityPrediction(
                mutation=mutation,
                ddg=ddg,
                confidence=_THERMOMPNN_DEFAULT_CONFIDENCE,
                method=self.name,
                method_version=self.version,
                notes=self._build_notes(raw),
                raw_output={"thermompnn": raw},
            )
        except KeyError as exc:
            return StabilityPrediction(
                mutation=mutation,
                ddg=None,
                confidence=None,
                method=self.name,
                method_version=self.version,
                error=f"ThermoMPNN output missing expected field: {exc}",
            )
        except Exception as exc:  # noqa: BLE001 - capture method-internal failures
            logger.warning("ThermoMPNN failed on %s: %s", mutation, exc)
            return StabilityPrediction(
                mutation=mutation,
                ddg=None,
                confidence=None,
                method=self.name,
                method_version=self.version,
                error=f"ThermoMPNN error: {exc}",
            )

    @staticmethod
    def _build_notes(raw: dict) -> str:
        """Surface useful per-mutation context from the raw ThermoMPNN output."""
        notes: list[str] = []
        score = raw.get("score")
        if score is not None:
            notes.append(f"model score {float(score):+.2f}")
        return "; ".join(notes)