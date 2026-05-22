"""FoldX method wrapper for stability prediction.

FoldX 5: Schymkowitz J et al. (2005) Nucleic Acids Res 33:W382-W388.
https://foldxsuite.crg.eu/

FoldX is an empirical force field that estimates folding free energy from
structure. ΔΔG is computed as Energy(mutant) - Energy(wild-type).

Workflow per mutation:
  1. RepairPDB on the input structure (one-time; cached per structure)
  2. BuildModel with the mutation specification
  3. Parse Dif_<pdb>.fxout for the ΔΔG value

Licensing: free for academic use, paid commercial license. Users must
register at https://foldxsuite.crg.eu/ and add the binary to PATH.

Ported and simplified from the foldx_wrapper.py in evosite3d
(github.com/romainstuder/evosite3d) by the same author.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from .base import Mutation, StabilityMethod, StabilityPrediction

logger = logging.getLogger(__name__)

# FoldX BuildModel output filename pattern:
#   Dif_<basename>.fxout, with a header row and one row per mutation set
_DIF_FILE_PATTERN = "Dif_{stem}.fxout"

# FoldX three-letter / one-letter amino acid mapping for error messages.
# FoldX itself uses one-letter codes in the mutation specification syntax.
_VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# FoldX confidence is roughly inversely proportional to local structure
# disorder. We don't expose a continuous confidence score yet; use a
# midrange default and let the consensus layer down-weight on disagreement.
_FOLDX_DEFAULT_CONFIDENCE = 0.75

# Cache RepairPDB results across mutations on the same structure.
# Keyed by (absolute_path, mtime) so we re-repair if the file changes.
_REPAIR_CACHE: dict[tuple[str, float], Path] = {}


class FoldX(StabilityMethod):
    """FoldX empirical force field stability predictor.

    Optional method. Requires manual installation of the FoldX binary
    (free for academic use). If unavailable, the skill falls back to
    RaSP and/or ThermoMPNN.
    """

    name = "foldx"

    def __init__(self, foldx_executable: str | None = None) -> None:
        """
        Args:
            foldx_executable: Path to the FoldX binary. If None, searches PATH
                for an executable named 'foldx', 'foldx5', or 'FoldX'.
        """
        super().__init__()
        self._executable: str | None = self._locate_executable(foldx_executable)
        self.version = self._probe_version() if self._executable else "not-installed"

    @staticmethod
    def _locate_executable(explicit_path: str | None) -> str | None:
        """Find a FoldX binary. Returns None if not found."""
        if explicit_path:
            return explicit_path if os.access(explicit_path, os.X_OK) else None
        for name in ("foldx", "foldx5", "FoldX"):
            found = shutil.which(name)
            if found:
                return found
        return None

    def _probe_version(self) -> str:
        """Run `foldx --version` (or equivalent) to capture the version string."""
        assert self._executable is not None
        try:
            # FoldX 5 prints version on --version; FoldX 4 prints on -h.
            # Try --version first, fall back to parsing -h output.
            result = subprocess.run(
                [self._executable, "--version"],
                capture_output=True, text=True, timeout=10,
            )
            output = (result.stdout + result.stderr).strip()
            match = re.search(r"FoldX\s+(\d+(?:\.\d+)*)", output)
            if match:
                return match.group(1)
            # Fall back to whatever was printed
            return output.split("\n")[0][:80] or "unknown"
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as exc:
            logger.warning("Could not probe FoldX version: %s", exc)
            return "unknown"

    def is_available(self) -> bool:
        return self._executable is not None

    def install_instructions(self) -> str:
        return (
            "FoldX is free for academic use but requires manual setup:\n"
            "\n"
            "  1. Register at https://foldxsuite.crg.eu/\n"
            "  2. Download the FoldX 5 binary for your OS (Linux, macOS, Windows)\n"
            "  3. Make it executable and add to PATH:\n"
            "       chmod +x foldx\n"
            "       export PATH=$PATH:/path/to/foldx5/\n"
            "  4. Verify:\n"
            "       foldx --version\n"
            "\n"
            "Commercial use requires a paid licence. See https://foldxsuite.crg.eu/\n"
            "\n"
            "Without FoldX, use --method rasp (default) or --method thermompnn,\n"
            "which are pip-installable and require no external setup.\n"
        )

    # ------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------

    def predict(
        self,
        structure_path: Path,
        mutations: list[Mutation],
    ) -> list[StabilityPrediction]:
        if not self.is_available():
            return [self._unavailable_prediction(m) for m in mutations]

        if not structure_path.exists():
            raise RuntimeError(f"Structure file not found: {structure_path}")

        # Repair the structure once; reuse across all mutations on this file.
        try:
            repaired_pdb = self._get_or_repair(structure_path)
        except RuntimeError as exc:
            return [
                StabilityPrediction(
                    mutation=m, ddg=None, confidence=None,
                    method=self.name, method_version=self.version,
                    error=f"FoldX RepairPDB failed: {exc}",
                )
                for m in mutations
            ]

        # Run BuildModel per mutation. Could be batched, but per-mutation
        # gives us cleaner per-result error handling.
        return [self._predict_one(repaired_pdb, m) for m in mutations]

    def _unavailable_prediction(self, mutation: Mutation) -> StabilityPrediction:
        return StabilityPrediction(
            mutation=mutation, ddg=None, confidence=None,
            method=self.name, method_version=self.version,
            error="FoldX binary not found. See install_instructions().",
        )

    def _predict_one(self, repaired_pdb: Path, mutation: Mutation) -> StabilityPrediction:
        """Run BuildModel for one mutation; return a StabilityPrediction."""
        validation_error = _validate_mutation(mutation)
        if validation_error:
            return StabilityPrediction(
                mutation=mutation, ddg=None, confidence=None,
                method=self.name, method_version=self.version,
                error=validation_error,
            )

        with tempfile.TemporaryDirectory(prefix="foldx_buildmodel_") as workdir:
            workdir_path = Path(workdir)
            try:
                self._copy_repaired_pdb(repaired_pdb, workdir_path)
                mutant_file = self._write_individual_list(workdir_path, mutation)
                self._run_build_model(workdir_path, repaired_pdb.name, mutant_file)
                ddg = self._parse_ddg(workdir_path, repaired_pdb.stem)
            except FoldXRuntimeError as exc:
                return StabilityPrediction(
                    mutation=mutation, ddg=None, confidence=None,
                    method=self.name, method_version=self.version,
                    error=str(exc),
                )

        return StabilityPrediction(
            mutation=mutation,
            ddg=ddg,
            confidence=_FOLDX_DEFAULT_CONFIDENCE,
            method=self.name,
            method_version=self.version,
            notes="BuildModel single-mutation run",
            raw_output={"ddg_kcal_mol": ddg},
        )

    # ------------------------------------------------------------------
    # Repair (cached per structure)
    # ------------------------------------------------------------------

    def _get_or_repair(self, structure_path: Path) -> Path:
        """Run RepairPDB once per structure; cache the result.

        Returns the path to the *_Repair.pdb file. The cached file lives in
        a stable temp location keyed by the input file's path and mtime.
        """
        abs_path = str(structure_path.resolve())
        mtime = structure_path.stat().st_mtime
        cache_key = (abs_path, mtime)
        if cache_key in _REPAIR_CACHE and _REPAIR_CACHE[cache_key].exists():
            return _REPAIR_CACHE[cache_key]

        repair_dir = Path(tempfile.mkdtemp(prefix="foldx_repair_"))
        local_copy = repair_dir / structure_path.name
        shutil.copy(structure_path, local_copy)
        self._run_repair(repair_dir, local_copy.name)

        repaired = repair_dir / f"{structure_path.stem}_Repair.pdb"
        if not repaired.exists():
            raise RuntimeError(f"RepairPDB completed but {repaired.name} not found")
        _REPAIR_CACHE[cache_key] = repaired
        return repaired

    def _run_repair(self, workdir: Path, pdb_filename: str) -> None:
        assert self._executable is not None
        cmd = [
            self._executable,
            "--command=RepairPDB",
            f"--pdb={pdb_filename}",
            f"--pdb-dir={workdir}",
            f"--output-dir={workdir}",
        ]
        self._run_foldx(cmd, label="RepairPDB")

    # ------------------------------------------------------------------
    # BuildModel
    # ------------------------------------------------------------------

    @staticmethod
    def _copy_repaired_pdb(repaired_pdb: Path, workdir: Path) -> None:
        shutil.copy(repaired_pdb, workdir / repaired_pdb.name)

    @staticmethod
    def _write_individual_list(workdir: Path, mutation: Mutation) -> str:
        """Create FoldX's individual_list.txt with the mutation specification.

        Format: one mutation set per line, ending with semicolon.
        Per-mutation syntax: <wt><chain><position><mt>
        Example: AA328S; means A in chain A at position 328 → S.
        """
        spec = f"{mutation.wt}{mutation.chain}{mutation.position}{mutation.mt};"
        list_file = workdir / "individual_list.txt"
        list_file.write_text(spec + "\n")
        return list_file.name

    def _run_build_model(self, workdir: Path, pdb_filename: str, mutant_file: str) -> None:
        assert self._executable is not None
        cmd = [
            self._executable,
            "--command=BuildModel",
            f"--pdb={pdb_filename}",
            f"--pdb-dir={workdir}",
            f"--mutant-file={mutant_file}",
            f"--output-dir={workdir}",
            "--numberOfRuns=1",
        ]
        self._run_foldx(cmd, label="BuildModel")

    def _run_foldx(self, cmd: list[str], *, label: str) -> None:
        """Execute a FoldX subprocess command, capturing failures with context."""
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)
        except subprocess.CalledProcessError as exc:
            raise FoldXRuntimeError(
                f"{label} failed (exit {exc.returncode}): {exc.stderr.strip()[:200]}"
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise FoldXRuntimeError(f"{label} timed out after {exc.timeout}s") from exc

    # ------------------------------------------------------------------
    # Output parsing
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_ddg(workdir: Path, pdb_stem: str) -> float:
        """Extract ΔΔG from Dif_<stem>.fxout.

        File format (FoldX 5):
            # comments...
            <header row with column names>
            <pdb_name>_1.pdb    <total ddg>    <breakdown columns>...

        The total ΔΔG is the second column of the first data row.
        """
        dif_file = workdir / _DIF_FILE_PATTERN.format(stem=pdb_stem)
        if not dif_file.exists():
            raise FoldXRuntimeError(f"Expected output {dif_file.name} not found")

        text = dif_file.read_text()
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Skip the header row (contains column names like "total energy")
            if line.lower().startswith("pdb") or "total energy" in line.lower():
                continue
            # First data row: <pdbname>\t<ddg>\t<other columns>
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                return float(parts[1])
            except ValueError:
                continue
        raise FoldXRuntimeError(f"Could not parse ΔΔG from {dif_file.name}")


# ----------------------------------------------------------------------------
# Validation helpers and exceptions
# ----------------------------------------------------------------------------

class FoldXRuntimeError(RuntimeError):
    """Raised internally when a FoldX subprocess or output parsing fails.

    Caught by _predict_one and converted into a per-mutation error.
    Not exported; users see StabilityPrediction.error instead.
    """


def _validate_mutation(mutation: Mutation) -> str | None:
    """Return an error string if the mutation is invalid, else None."""
    if mutation.mt == "-":
        return "Deletions not supported in v1"
    if mutation.wt not in _VALID_AA:
        return f"Invalid wild-type amino acid: {mutation.wt!r}"
    if mutation.mt not in _VALID_AA:
        return f"Invalid mutant amino acid: {mutation.mt!r}"
    if mutation.wt == mutation.mt:
        return "Wild-type and mutant residues are identical; nothing to predict"
    if not mutation.chain or len(mutation.chain) != 1:
        return f"Chain must be a single character; got {mutation.chain!r}"
    if mutation.position < 1:
        return f"Position must be >= 1; got {mutation.position}"
    return None
