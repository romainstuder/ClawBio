"""Input/output utilities for stability-predictor.

Responsibilities:
  - Parse and validate the mutations JSON input
  - Inspect structure files (PDB/CIF) for chain/residue info
  - Write result.json, predictions.json, and the reproducibility bundle
  - Compute SHA-256 checksums of inputs and outputs

Kept dependency-light: structure inspection uses Biopython if available but
degrades gracefully (returns a placeholder chain summary) if not. Biopython
is listed as a hard dep in SKILL.md, but unit tests shouldn't require it.
"""

from __future__ import annotations

import hashlib
import json
import logging
import platform
import shlex
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .consensus import ConsensusResult
from .methods.base import Mutation

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------------
# Mutation input
# ----------------------------------------------------------------------------

class MutationInputError(ValueError):
    """Raised when mutation JSON is malformed or invalid."""


def load_mutations(path: Path) -> list[Mutation]:
    """Load and validate mutations from JSON.

    Expected format:
        {
            "chain": "A",
            "mutations": [
                {"position": 328, "wt": "A", "mt": "S"},
                {"position": 508, "wt": "F", "mt": "A", "chain": "B"}
            ]
        }

    Top-level "chain" is the default; each mutation may override it.

    Returns:
        List of validated Mutation objects in input order.

    Raises:
        MutationInputError: If JSON is malformed or mutations are invalid.
    """
    if not path.exists():
        raise MutationInputError(f"Mutations file not found: {path}")

    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise MutationInputError(f"Invalid JSON in {path.name}: {exc}") from exc

    if not isinstance(data, dict):
        raise MutationInputError(f"Top-level JSON must be an object, got {type(data).__name__}")

    default_chain = data.get("chain")
    raw_mutations = data.get("mutations")
    if not isinstance(raw_mutations, list):
        raise MutationInputError("'mutations' must be a list")
    if not raw_mutations:
        raise MutationInputError("'mutations' list is empty")

    mutations: list[Mutation] = []
    for i, raw in enumerate(raw_mutations):
        if not isinstance(raw, dict):
            raise MutationInputError(f"Mutation #{i} must be an object, got {type(raw).__name__}")
        try:
            mutations.append(_parse_one_mutation(raw, default_chain, i))
        except MutationInputError:
            raise
        except (TypeError, ValueError) as exc:
            raise MutationInputError(f"Mutation #{i} invalid: {exc}") from exc
    return mutations


def _parse_one_mutation(raw: dict, default_chain: str | None, index: int) -> Mutation:
    chain = raw.get("chain", default_chain)
    if not chain:
        raise MutationInputError(
            f"Mutation #{index}: no chain specified and no default 'chain' at top level"
        )

    for required in ("position", "wt", "mt"):
        if required not in raw:
            raise MutationInputError(f"Mutation #{index} missing required field: {required!r}")

    return Mutation(
        chain=str(chain),
        position=int(raw["position"]),
        wt=str(raw["wt"]).upper(),
        mt=str(raw["mt"]).upper(),
    )


# ----------------------------------------------------------------------------
# Structure inspection
# ----------------------------------------------------------------------------

@dataclass
class StructureSummary:
    """Lightweight summary of an input structure for the report header."""

    chains: list[str]
    residue_count: int
    method_used: str  # "biopython" or "fallback"

    def human_readable(self, focus_chain: str | None = None) -> str:
        chain_str = ", ".join(self.chains) if self.chains else "unknown chains"
        if focus_chain and focus_chain in self.chains:
            chain_str = f"chain {focus_chain} (of {chain_str})"
        elif len(self.chains) == 1:
            chain_str = f"chain {self.chains[0]}"
        return f"{chain_str}, {self.residue_count} residues"


def inspect_structure(structure_path: Path) -> StructureSummary:
    """Return a lightweight summary of the structure for reporting purposes.

    Uses Biopython if available; otherwise returns a fallback summary based
    on parsing the file directly for ATOM/HETATM records (PDB only).
    """
    try:
        return _inspect_with_biopython(structure_path)
    except ImportError:
        logger.info("Biopython unavailable; using fallback PDB parser")
        return _inspect_fallback(structure_path)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Biopython structure inspection failed: %s. Falling back.", exc)
        return _inspect_fallback(structure_path)


def _inspect_with_biopython(structure_path: Path) -> StructureSummary:
    """Use Biopython's structure parser."""
    from Bio.PDB import MMCIFParser, PDBParser  # type: ignore

    suffix = structure_path.suffix.lower()
    parser: Any
    if suffix in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("input", str(structure_path))
    chains: list[str] = []
    residue_count = 0
    for model in structure:
        for chain in model:
            chains.append(chain.id)
            residue_count += sum(1 for _ in chain.get_residues())
        break  # first model only
    return StructureSummary(chains=chains, residue_count=residue_count, method_used="biopython")


def _inspect_fallback(structure_path: Path) -> StructureSummary:
    """Crude PDB-only fallback: count unique chain IDs and CA atoms.

    Not accurate for CIF; if the file isn't a PDB, returns empty summary.
    """
    if structure_path.suffix.lower() not in (".pdb", ".ent"):
        return StructureSummary(chains=[], residue_count=0, method_used="fallback-empty")

    chains: set[str] = set()
    ca_count = 0
    with structure_path.open() as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain_id = line[21:22].strip()
            atom_name = line[12:16].strip()
            if chain_id:
                chains.add(chain_id)
            if atom_name == "CA":
                ca_count += 1
    return StructureSummary(chains=sorted(chains), residue_count=ca_count, method_used="fallback")


# ----------------------------------------------------------------------------
# Output writing
# ----------------------------------------------------------------------------

def write_result_json(
    output_dir: Path,
    *,
    consensus_results: list[ConsensusResult],
    structure_path: Path,
    methods_requested: list[str],
    methods_available: list[str],
) -> None:
    """Write result.json — the ClawBio-standard top-level summary."""
    n = len(consensus_results)
    n_destabilizing = sum(
        1 for c in consensus_results if c.direction.value == "destabilizing"
    )
    n_neutral = sum(1 for c in consensus_results if c.direction.value == "neutral")
    n_stabilizing = sum(
        1 for c in consensus_results if c.direction.value == "stabilizing"
    )
    n_flagged = sum(1 for c in consensus_results if c.flags)

    payload = {
        "skill": "stability-predictor",
        "version": "0.1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "input": {
            "structure": str(structure_path.name),
            "n_mutations": n,
        },
        "methods": {
            "requested": methods_requested,
            "available": methods_available,
            "skipped": sorted(set(methods_requested) - set(methods_available)),
        },
        "summary": {
            "n_destabilizing": n_destabilizing,
            "n_neutral": n_neutral,
            "n_stabilizing": n_stabilizing,
            "n_flagged_for_review": n_flagged,
        },
    }
    _write_json(output_dir / "result.json", payload)


def write_predictions_json(
    output_dir: Path,
    consensus_results: list[ConsensusResult],
) -> None:
    """Write predictions.json — per-mutation, per-method machine-readable details."""
    payload = {"predictions": [c.to_dict() for c in consensus_results]}
    _write_json(output_dir / "predictions.json", payload)


def write_reproducibility_bundle(
    output_dir: Path,
    *,
    argv: list[str],
    methods_available: list[str],
    method_versions: dict[str, str],
    input_paths: list[Path],
) -> None:
    """Write the reproducibility/ subdirectory.

    Contents:
      commands.sh         — exact replay command
      environment.yml     — Python + method versions
      input_checksum.txt  — SHA-256 of inputs
      output_checksum.txt — SHA-256 of outputs (written by caller after this)
    """
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    # commands.sh
    cmd_line = " ".join(shlex.quote(a) for a in argv)
    commands_sh = (
        "#!/usr/bin/env bash\n"
        "# Auto-generated replay script for stability-predictor.\n"
        f"# Generated: {datetime.now(timezone.utc).isoformat()}\n"
        "set -euo pipefail\n\n"
        f"{cmd_line}\n"
    )
    (repro_dir / "commands.sh").write_text(commands_sh)
    (repro_dir / "commands.sh").chmod(0o755)

    # environment.yml (a minimal YAML; full env capture is overkill for v1)
    env_lines = [
        "# Stability-predictor environment snapshot",
        f"python: {sys.version.split()[0]}",
        f"platform: {platform.platform()}",
        "methods:",
    ]
    for method_name in methods_available:
        version = method_versions.get(method_name, "unknown")
        env_lines.append(f"  {method_name}: {version}")
    (repro_dir / "environment.yml").write_text("\n".join(env_lines) + "\n")

    # input_checksum.txt
    checksum_lines = []
    for path in input_paths:
        if path.exists():
            checksum_lines.append(f"{sha256_file(path)}  {path.name}")
    (repro_dir / "input_checksum.txt").write_text("\n".join(checksum_lines) + "\n")


def write_output_checksums(output_dir: Path) -> None:
    """Compute and write SHA-256 for all files in output_dir except checksums themselves.

    Call this *after* all other outputs have been written.
    """
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    checksum_file = repro_dir / "output_checksum.txt"

    lines = []
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file():
            continue
        if path == checksum_file:
            continue
        relative = path.relative_to(output_dir)
        lines.append(f"{sha256_file(path)}  {relative}")
    checksum_file.write_text("\n".join(lines) + "\n")


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n")


def sha256_file(path: Path, chunk_size: int = 65536) -> str:
    """SHA-256 hex digest of a file's contents."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()
