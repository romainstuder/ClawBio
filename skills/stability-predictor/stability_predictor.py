#!/usr/bin/env python3
"""Stability Predictor — predict ΔΔG of point mutations on a protein structure.

ClawBio skill. Wraps RaSP (default), ThermoMPNN, and FoldX (optional).
See SKILL.md for full specification, README.md for quick start.

Usage:
    # Default: RaSP on user data
    python stability_predictor.py --structure protein.pdb --mutations mutations.json --output results/

    # Specific method
    python stability_predictor.py --structure ... --mutations ... --method foldx --output ...

    # All available methods + consensus
    python stability_predictor.py --structure ... --mutations ... --method all --output ...

    # Clinical demo (CFTR F508del)
    python stability_predictor.py --demo --output /tmp/sp_demo

    # Scientific demo (RubisCO, validates against Studer et al. 2014)
    python stability_predictor.py --demo --demo-set rubisco --output /tmp/sp_demo
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from stability_predictor_core.consensus import aggregate
from stability_predictor_core.io import (
    MutationInputError,
    inspect_structure,
    load_mutations,
    write_output_checksums,
    write_predictions_json,
    write_reproducibility_bundle,
    write_result_json,
)
from stability_predictor_core.methods.base import Mutation, StabilityMethod
from stability_predictor_core.methods.foldx import FoldX
from stability_predictor_core.methods.precomputed import Precomputed
from stability_predictor_core.methods.rasp import RaSP
from stability_predictor_core.plotting import render_figures
from stability_predictor_core.report import render_report

# ThermoMPNN import is best-effort; the skill should still run if its
# module isn't present (e.g., during early development).
try:
    from stability_predictor_core.methods.thermompnn import ThermoMPNN  # type: ignore
    _THERMOMPNN_AVAILABLE = True
except ImportError:
    ThermoMPNN = None  # type: ignore
    _THERMOMPNN_AVAILABLE = False

logger = logging.getLogger("stability_predictor")

# Method registry. Add new methods here and they become available via --method.
_METHOD_REGISTRY: dict[str, type[StabilityMethod]] = {"rasp": RaSP, "foldx": FoldX}
if _THERMOMPNN_AVAILABLE:
    _METHOD_REGISTRY["thermompnn"] = ThermoMPNN  # type: ignore

DEFAULT_METHOD = "rasp"
METHOD_ALL = "all"

# Demo data lives next to this script under demo_data/.
_SCRIPT_DIR = Path(__file__).resolve().parent
_DEMO_DIR = _SCRIPT_DIR / "demo_data"
_DEMO_SETS = {
    "t4lysozyme": {
        "structure": _DEMO_DIR / "t4lysozyme.pdb",
        "mutations": _DEMO_DIR / "t4lysozyme_mutations.json",
        "caches": {
            "rasp": _DEMO_DIR / "t4lysozyme_predictions_rasp.json",
            "thermompnn": _DEMO_DIR / "t4lysozyme_predictions_thermompnn.json",
            "foldx": _DEMO_DIR / "t4lysozyme_predictions_foldx.json",
        },
        "description": (
            "T4 lysozyme (PDB 2LZM): L99A cavity-creating mutation "
            "(Eriksson et al. 1992) plus T157I and T26S benchmarks from "
            "the Matthews lab dataset."
        ),
    },
    "p53": {
        "structure": _DEMO_DIR / "p53.pdb",
        "mutations": _DEMO_DIR / "p53_mutations.json",
        "caches": {
            "rasp": _DEMO_DIR / "p53_predictions_rasp.json",
            "thermompnn": _DEMO_DIR / "p53_predictions_thermompnn.json",
            "foldx": _DEMO_DIR / "p53_predictions_foldx.json",
        },
        "description": (
            "p53 DNA-binding domain (PDB 2XWR): Y220C destabilising "
            "cancer mutation (Joerger & Fersht 2007)."
        ),
    },
}

EXIT_OK = 0
EXIT_USAGE = 2
EXIT_INPUT_ERROR = 3
EXIT_METHOD_UNAVAILABLE = 4
EXIT_RUNTIME_ERROR = 5


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="stability-predictor",
        description="Predict ΔΔG of point mutations using RaSP, ThermoMPNN, and/or FoldX.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Methods:\n"
            "  rasp        Open source, pip-installable, fast (default)\n"
            "  thermompnn  Open source, pip-installable, SOTA on Megascale dataset\n"
            "  foldx       Optional; free for academics; requires manual install\n"
            "  all         Run every available method and compute consensus\n"
            "\n"
            "Demo sets:\n"
            "  t4lysozyme  L99A / T157I / T26S — Matthews-lab benchmarks (default)\n"
            "  p53         Y220C destabilising cancer mutation (Joerger & Fersht 2007)\n"
        ),
    )
    parser.add_argument("--structure", type=Path, help="PDB or CIF structure file")
    parser.add_argument("--mutations", type=Path, help="JSON file listing mutations")
    parser.add_argument(
        "--method",
        default=DEFAULT_METHOD,
        choices=sorted(_METHOD_REGISTRY.keys()) + [METHOD_ALL],
        help="Prediction method (default: rasp)",
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="Output directory (will be created)"
    )
    parser.add_argument(
        "--demo", action="store_true", help="Run on built-in demo data (skips --structure/--mutations)"
    )
    parser.add_argument(
        "--demo-set",
        default="t4lysozyme",
        choices=sorted(_DEMO_SETS.keys()),
        help="Which demo dataset to use with --demo (default: t4lysozyme)",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Verbose logging"
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    # Resolve inputs (demo or user-supplied)
    try:
        structure_path, mutations_path = _resolve_inputs(args)
    except _UsageError as exc:
        logger.error(str(exc))
        return EXIT_USAGE

    # Load and validate mutations
    try:
        mutations = load_mutations(mutations_path)
    except MutationInputError as exc:
        logger.error("Mutation input error: %s", exc)
        return EXIT_INPUT_ERROR

    if not structure_path.exists():
        logger.error("Structure file not found: %s", structure_path)
        return EXIT_INPUT_ERROR

    logger.info(
        "Loaded %d mutation(s) on %s", len(mutations), structure_path.name
    )

    # Resolve which methods to run. In demo mode we substitute cached
    # predictions so the run completes in seconds without external models.
    methods_requested = _resolve_method_selection(args.method)
    if args.demo:
        methods, missing = _instantiate_demo_methods(
            methods_requested, _DEMO_SETS[args.demo_set]["caches"]
        )
    else:
        methods, missing = _instantiate_methods(methods_requested)
    if not methods:
        logger.error(
            "None of the requested methods are available: %s", methods_requested
        )
        for name in missing:
            logger.error("\n--- %s install instructions ---\n%s",
                         name, _METHOD_REGISTRY[name]().install_instructions())
        return EXIT_METHOD_UNAVAILABLE

    methods_available = [m.name for m in methods]
    if missing:
        logger.warning(
            "Skipping unavailable methods: %s (run with --verbose to see install instructions)",
            ", ".join(missing),
        )
        if args.verbose:
            for name in missing:
                logger.info("\n--- %s ---\n%s",
                            name, _METHOD_REGISTRY[name]().install_instructions())

    # Prepare output directory
    output_dir: Path = args.output
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)

    # Run predictions
    try:
        predictions_per_mutation = _run_all_methods(methods, structure_path, mutations)
    except RuntimeError as exc:
        logger.error("Prediction run failed: %s", exc)
        return EXIT_RUNTIME_ERROR

    # Aggregate into consensus
    consensus_results = aggregate(predictions_per_mutation)

    # Inspect structure for the report header
    structure_summary = inspect_structure(structure_path)
    focus_chain = mutations[0].chain if mutations else None
    chain_summary = structure_summary.human_readable(focus_chain=focus_chain)

    # Write outputs
    _write_outputs(
        output_dir=output_dir,
        consensus_results=consensus_results,
        structure_path=structure_path,
        mutations_path=mutations_path,
        chain_summary=chain_summary,
        methods_requested=methods_requested,
        methods_available=methods_available,
        method_versions={m.name: m.version for m in methods},
        method=args.method,
    )

    # Console summary
    _print_summary(consensus_results, output_dir)
    return EXIT_OK


# ----------------------------------------------------------------------------
# Step helpers
# ----------------------------------------------------------------------------

class _UsageError(Exception):
    """Internal: bad argument combination."""


def _resolve_inputs(args: argparse.Namespace) -> tuple[Path, Path]:
    """Return (structure_path, mutations_path), resolving demo flags if set."""
    if args.demo:
        demo = _DEMO_SETS[args.demo_set]
        logger.info("Demo mode: %s — %s", args.demo_set, demo["description"])
        return demo["structure"], demo["mutations"]
    if not args.structure or not args.mutations:
        raise _UsageError(
            "Either --demo or both --structure and --mutations are required"
        )
    return args.structure, args.mutations


def _resolve_method_selection(requested: str) -> list[str]:
    """Expand 'all' into the full list of registered methods."""
    if requested == METHOD_ALL:
        return sorted(_METHOD_REGISTRY.keys())
    return [requested]


def _instantiate_demo_methods(
    requested: list[str],
    caches: dict[str, Path],
) -> tuple[list[StabilityMethod], list[str]]:
    """In demo mode, back each requested method with a Precomputed cache.

    Missing caches are reported the same way unavailable real methods are.
    """
    available: list[StabilityMethod] = []
    missing: list[str] = []
    for name in requested:
        cache_path = caches.get(name)
        if cache_path is None:
            missing.append(name)
            continue
        method = Precomputed(cache_path, label=name)
        if method.is_available():
            available.append(method)
        else:
            missing.append(name)
    return available, missing


def _instantiate_methods(
    requested: list[str],
) -> tuple[list[StabilityMethod], list[str]]:
    """Instantiate each requested method; partition into available and missing."""
    available: list[StabilityMethod] = []
    missing: list[str] = []
    for name in requested:
        method_class = _METHOD_REGISTRY[name]
        method = method_class()
        if method.is_available():
            available.append(method)
        else:
            missing.append(name)
    return available, missing


def _run_all_methods(
    methods: list[StabilityMethod],
    structure_path: Path,
    mutations: list[Mutation],
) -> list[list]:
    """Run each method on all mutations; reshape into per-mutation lists.

    Returns:
        Outer list: one entry per mutation, in input order.
        Inner list: one StabilityPrediction per method.
    """
    per_method: dict[str, list] = {}
    for method in methods:
        logger.info("Running %s (v%s) on %d mutation(s)...",
                    method.name, method.version, len(mutations))
        per_method[method.name] = method.predict(structure_path, mutations)

    # Reshape: list-of-lists keyed by mutation index
    per_mutation: list[list] = []
    for i in range(len(mutations)):
        per_mutation.append([per_method[m.name][i] for m in methods])
    return per_mutation


def _write_outputs(
    *,
    output_dir: Path,
    consensus_results: list,
    structure_path: Path,
    mutations_path: Path,
    chain_summary: str,
    methods_requested: list[str],
    methods_available: list[str],
    method_versions: dict[str, str],
    method: str,
) -> None:
    """Write all output files: figures, report.md, JSONs, reproducibility bundle, checksums."""
    # Figures (best-effort: silently skipped if matplotlib is unavailable)
    render_figures(consensus_results, methods_available, output_dir)

    # Report
    md = render_report(
        consensus_results=consensus_results,
        structure_path=structure_path,
        chain_summary=chain_summary,
        methods_requested=methods_requested,
        methods_available=methods_available,
        output_dir=output_dir,
    )
    (output_dir / "report.md").write_text(md)

    # Machine-readable
    write_result_json(
        output_dir,
        consensus_results=consensus_results,
        structure_path=structure_path,
        methods_requested=methods_requested,
        methods_available=methods_available,
    )
    write_predictions_json(output_dir, consensus_results)

    # Reproducibility: rebuild a portable replay invocation from parsed args
    # rather than echoing sys.argv (which depends on how the CLI was launched).
    replay_argv = [
        "python", "stability_predictor.py",
        "--structure", str(structure_path),
        "--mutations", str(mutations_path),
        "--method", method,
        "--output", str(output_dir),
    ]
    write_reproducibility_bundle(
        output_dir,
        argv=replay_argv,
        methods_available=methods_available,
        method_versions=method_versions,
        input_paths=[structure_path, mutations_path],
    )
    # Must come last: hashes every file in the output dir
    write_output_checksums(output_dir)


def _print_summary(consensus_results: list, output_dir: Path) -> None:
    """One-screen summary printed to stdout after a successful run."""
    n = len(consensus_results)
    n_destab = sum(1 for c in consensus_results if c.direction.value == "destabilizing")
    n_neutral = sum(1 for c in consensus_results if c.direction.value == "neutral")
    n_stab = sum(1 for c in consensus_results if c.direction.value == "stabilizing")
    n_unknown = sum(1 for c in consensus_results if c.direction.value == "unknown")
    n_flagged = sum(1 for c in consensus_results if c.flags)

    print()
    print(f"Predicted {n} mutation(s):")
    print(f"  destabilizing: {n_destab}")
    print(f"  neutral:       {n_neutral}")
    print(f"  stabilizing:   {n_stab}")
    if n_unknown:
        print(f"  unknown:       {n_unknown}  (all methods failed)")
    if n_flagged:
        print(f"  flagged for review: {n_flagged}")
    print()
    print(f"Outputs written to: {output_dir}")
    print(f"  report.md          — human-readable summary")
    print(f"  result.json        — top-level machine-readable")
    print(f"  predictions.json   — per-mutation per-method detail")
    print(f"  reproducibility/   — replay bundle + checksums")
    print()


if __name__ == "__main__":
    sys.exit(main())
