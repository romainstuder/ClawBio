"""Precomputed-cache stability method (for fast demos and offline tests).

Returns ΔΔG values from a JSON lookup table keyed by (chain, position, wt, mt).
Used for `--demo` so users don't need RaSP / ThermoMPNN / FoldX installed to
see the skill work end-to-end.

The cached values are literature-grounded estimates; they are not new
experimental measurements and should not be used as a benchmark. The intent
is purely to exercise the pipeline and produce a representative report.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

from .base import Mutation, StabilityMethod, StabilityPrediction

logger = logging.getLogger(__name__)


class Precomputed(StabilityMethod):
    """Looks up ΔΔG from a per-method JSON cache.

    The cache is a JSON file with shape::

        {
          "method": "rasp",
          "method_version": "1.2.0 (cached)",
          "predictions": [
            {"chain": "A", "position": 508, "wt": "F", "mt": "A",
             "ddg": 3.21, "notes": "destabilizing"},
            ...
          ]
        }

    Looking up a mutation not in the cache returns an explicit error
    (the prediction "succeeded" property is False).
    """

    name = "precomputed"

    def __init__(self, cache_path: Path | str, label: str | None = None) -> None:
        super().__init__()
        self._cache_path = Path(cache_path)
        self._label = label
        self._index: dict[tuple[str, int, str, str], dict] = {}
        self._cache_method_name = "precomputed"
        self._cache_method_version = "cached-0.1"
        self._load_error: str | None = None
        self._load_cache()

    # ------------------------------------------------------------------
    # Cache loading
    # ------------------------------------------------------------------

    def _load_cache(self) -> None:
        if not self._cache_path.exists():
            self._load_error = f"Cache file not found: {self._cache_path}"
            return
        try:
            data = json.loads(self._cache_path.read_text())
        except json.JSONDecodeError as exc:
            self._load_error = f"Invalid JSON in cache: {exc}"
            return

        self._cache_method_name = data.get("method", self._cache_method_name)
        self._cache_method_version = data.get("method_version", self._cache_method_version)
        predictions = data.get("predictions", [])
        if not isinstance(predictions, list):
            self._load_error = "'predictions' must be a list"
            return

        for entry in predictions:
            try:
                key = (
                    str(entry["chain"]),
                    int(entry["position"]),
                    str(entry["wt"]).upper(),
                    str(entry["mt"]).upper(),
                )
            except (KeyError, TypeError, ValueError) as exc:
                logger.warning("Skipping malformed cache entry %r: %s", entry, exc)
                continue
            self._index[key] = entry

        # Public name reflects what the cache claims to represent (e.g. "rasp"),
        # so reports and consensus tables read naturally.
        self.name = self._label or self._cache_method_name
        self.version = self._cache_method_version

    # ------------------------------------------------------------------
    # StabilityMethod interface
    # ------------------------------------------------------------------

    def is_available(self) -> bool:
        return self._load_error is None and bool(self._index)

    def install_instructions(self) -> str:
        if self._load_error:
            return (
                f"Precomputed cache unavailable: {self._load_error}\n"
                "This method ships with demo data and is normally only used "
                "via --demo. Check that demo_data/*_predictions_*.json exists "
                "and is readable.\n"
            )
        return "Precomputed cache is bundled with the skill; no setup needed.\n"

    def predict(
        self,
        structure_path: Path,
        mutations: list[Mutation],
    ) -> list[StabilityPrediction]:
        if not self.is_available():
            return [self._missing_prediction(m, self._load_error or "cache unavailable")
                    for m in mutations]

        results: list[StabilityPrediction] = []
        for mutation in mutations:
            entry = self._index.get(
                (mutation.chain, mutation.position, mutation.wt, mutation.mt)
            )
            if entry is None:
                results.append(self._missing_prediction(
                    mutation,
                    f"no cached prediction for {mutation}; "
                    f"regenerate the cache from demo_data/",
                ))
                continue
            results.append(StabilityPrediction(
                mutation=mutation,
                ddg=float(entry["ddg"]),
                confidence=float(entry.get("confidence", 0.7)),
                method=self.name,
                method_version=self.version,
                notes=str(entry.get("notes", "from cached predictions")),
                raw_output={"cache_path": str(self._cache_path)},
            ))
        return results

    def _missing_prediction(self, mutation: Mutation, error: str) -> StabilityPrediction:
        return StabilityPrediction(
            mutation=mutation, ddg=None, confidence=None,
            method=self.name, method_version=self.version,
            error=error,
        )