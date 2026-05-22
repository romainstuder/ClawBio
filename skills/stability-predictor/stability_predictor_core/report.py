"""Markdown report generation for stability predictions.

Produces `report.md` — the primary user-facing deliverable. Structure mirrors
the spec in SKILL.md and the example in README.md.

Design rules:
  - No HTML; pure Markdown for max portability
  - No emoji-as-data (only as section decoration, sparingly)
  - Tables are the unit of communication
  - Every claim cites a method + version
  - Limitations and reproducibility sections are mandatory, not optional
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from .consensus import ConsensusResult, Direction, Confidence

# Citations are defined centrally so we can keep them consistent across reports
# and update them in one place if a method's preferred citation changes.
_CITATIONS = {
    "rasp": "Blaabjerg LM et al. (2023) eLife 12:e82593. doi:10.7554/eLife.82593",
    "thermompnn": "Dieckhaus H et al. (2024) PNAS 121(6):e2314853121. doi:10.1073/pnas.2314853121",
    "foldx": "Schymkowitz J et al. (2005) Nucleic Acids Res 33(W):W382-W388. doi:10.1093/nar/gki387",
}

_VALIDATION_CITATION = (
    "Studer RA, Christin P-A, Williams MA, Orengo CA (2014) "
    "Stability-activity tradeoffs constrain the adaptive evolution of RubisCO. "
    "PNAS 111(6):2223-2228. doi:10.1073/pnas.1310811111"
)


def render_report(
    *,
    consensus_results: list[ConsensusResult],
    structure_path: Path,
    chain_summary: str,
    methods_requested: list[str],
    methods_available: list[str],
    output_dir: Path,
) -> str:
    """Render the full Markdown report.

    Args:
        consensus_results: One per mutation.
        structure_path: Input structure path (used in header).
        chain_summary: Human-readable chain info (e.g., "chain A, 1480 residues").
        methods_requested: Methods the user asked for.
        methods_available: Methods that were actually runnable on this machine.
        output_dir: Where the report will live; used for relative figure paths.

    Returns:
        Markdown report as a string. Caller writes to disk.
    """
    sections: list[str] = []
    sections.append(_render_header(structure_path, chain_summary, methods_requested,
                                    methods_available, len(consensus_results)))
    sections.append(_render_summary_table(consensus_results, methods_available))
    sections.append(_render_interpretation_guide())
    sections.append(_render_figures(output_dir))
    sections.append(_render_method_agreement(consensus_results, output_dir))
    sections.append(_render_flagged_mutations(consensus_results))
    sections.append(_render_limitations(consensus_results))
    sections.append(_render_methods_section(methods_available))
    sections.append(_render_references(methods_available))
    sections.append(_render_reproducibility(output_dir))
    sections.append(_render_disclaimer())
    return "\n\n".join(sections) + "\n"


# ----------------------------------------------------------------------------
# Section renderers (each returns a self-contained Markdown block)
# ----------------------------------------------------------------------------

def _render_header(
    structure_path: Path,
    chain_summary: str,
    methods_requested: list[str],
    methods_available: list[str],
    n_mutations: int,
) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    skipped = sorted(set(methods_requested) - set(methods_available))
    lines = [
        "# Stability Prediction Report",
        "",
        "## Input",
        "",
        f"- **Structure**: `{structure_path.name}` ({chain_summary})",
        f"- **Mutations**: {n_mutations}",
        f"- **Methods requested**: {', '.join(methods_requested) or 'none'}",
        f"- **Methods used**: {', '.join(methods_available) or 'none'}",
    ]
    if skipped:
        lines.append(
            f"- **Methods skipped**: {', '.join(skipped)} (not available on this system; "
            "see Methods section for install instructions)"
        )
    lines.append(f"- **Generated**: {timestamp}")
    return "\n".join(lines)


def _render_summary_table(
    consensus_results: list[ConsensusResult],
    methods_available: list[str],
) -> str:
    if not consensus_results:
        return "## Summary\n\nNo mutations predicted."

    header_cells = ["Mutation"] + [f"{m} ΔΔG" for m in methods_available] + [
        "Consensus", "Direction", "Confidence"
    ]
    rows = ["| " + " | ".join(header_cells) + " |",
            "|" + "|".join(["---"] * len(header_cells)) + "|"]

    for cr in consensus_results:
        per_method_values = []
        method_results = {p.method: p for p in cr.per_method}
        for method_name in methods_available:
            pred = method_results.get(method_name)
            if pred is None or not pred.succeeded:
                per_method_values.append("—")
            else:
                per_method_values.append(f"{pred.ddg:+.2f}")  # type: ignore[union-attr]

        consensus_cell = (
            f"{cr.consensus_ddg:+.2f}" if cr.consensus_ddg is not None else "—"
        )
        direction_cell = _format_direction(cr.direction)
        confidence_cell = _format_confidence(cr)

        row_cells = [str(cr.mutation)] + per_method_values + [
            consensus_cell, direction_cell, confidence_cell
        ]
        rows.append("| " + " | ".join(row_cells) + " |")

    return "## Summary\n\nΔΔG values in kcal/mol. Positive = destabilizing.\n\n" + "\n".join(rows)


def _render_interpretation_guide() -> str:
    return (
        "## Interpretation Guide\n\n"
        "| ΔΔG (kcal/mol) | Meaning |\n"
        "|---|---|\n"
        "| > +1.0 | **Destabilizing** — protein likely misfolds or has reduced fold stability |\n"
        "| -1.0 to +1.0 | **Neutral** — no major effect on folding stability |\n"
        "| < -1.0 | **Stabilizing** — improves folding stability |\n\n"
        "These thresholds approximate the energy of one hydrogen bond and follow "
        "common conventions in the directed-evolution literature."
    )


def _render_figures(output_dir: Path) -> str:
    """Embed any of the standard figures that were generated this run."""
    figures = [
        ("ddg_heatmap.png", "ΔΔG heatmap (mutations × methods)"),
        ("per_mutation_bars.png", "Per-mutation ΔΔG bars"),
    ]
    present = [
        (name, caption)
        for name, caption in figures
        if (output_dir / "figures" / name).exists()
    ]
    if not present:
        return "## Figures\n\nNo figures generated (matplotlib unavailable, or insufficient data)."
    lines = ["## Figures", ""]
    for name, caption in present:
        rel = Path("figures") / name
        lines.append(f"![{caption}]({rel})")
        lines.append("")
        lines.append(f"*{caption}*")
        lines.append("")
    return "\n".join(lines).rstrip()


def _render_method_agreement(
    consensus_results: list[ConsensusResult],
    output_dir: Path,
) -> str:
    n_total = len(consensus_results)
    if n_total == 0:
        return ""

    n_agree_direction = sum(1 for c in consensus_results if c.methods_agreed_direction)
    n_agree_magnitude = sum(1 for c in consensus_results if c.methods_agreed_magnitude)
    n_high_conf = sum(1 for c in consensus_results if c.confidence == Confidence.HIGH)

    figure_rel = Path("figures") / "method_agreement.png"
    figure_block = (
        f"![Method agreement scatter]({figure_rel})\n\n"
        if (output_dir / figure_rel).exists()
        else ""
    )

    return (
        "## Method Agreement\n\n"
        f"- **Direction agreement**: {n_agree_direction}/{n_total} mutations\n"
        f"- **Magnitude agreement** (within 1.0 kcal/mol): {n_agree_magnitude}/{n_total}\n"
        f"- **High-confidence calls**: {n_high_conf}/{n_total}\n\n"
        f"{figure_block}"
        "Disagreement between methods is informative, not a failure. Physics-based "
        "(FoldX) and ML-based (RaSP, ThermoMPNN) methods make different assumptions; "
        "where they diverge, the prediction warrants closer inspection."
    )


def _render_flagged_mutations(consensus_results: list[ConsensusResult]) -> str:
    flagged = [c for c in consensus_results if c.flags]
    if not flagged:
        return "## Flags\n\nNo mutations flagged for review."

    lines = ["## Flags\n\nThe following mutations have notes worth reviewing:"]
    for cr in flagged:
        lines.append(f"\n### {cr.mutation}")
        for flag in cr.flags:
            lines.append(f"- {flag}")
    return "\n".join(lines)


def _render_limitations(consensus_results: list[ConsensusResult]) -> str:
    return (
        "## Limitations\n\n"
        "- Predictions assume a **rigid backbone**. Mutations affecting hinge motions, "
        "loop flexibility, or large conformational rearrangements may be under-predicted.\n"
        "- All three methods predict **folding stability only**. Activity, binding affinity, "
        "and allosteric effects are not addressed.\n"
        "- ML methods (RaSP, ThermoMPNN) reflect their training distributions; predictions "
        "on protein families heavily under-represented in PDB may be less reliable.\n"
        "- FoldX is calibrated against directed-evolution datasets dominated by "
        "small, well-folded proteins.\n"
        "- Deletions, insertions, and multi-residue changes are **not supported** in v1.\n"
        "- Predictions are hypotheses, not experimental measurements."
    )


def _render_methods_section(methods_available: list[str]) -> str:
    descriptions = {
        "rasp": (
            "**RaSP** — CNN trained on Rosetta ΔΔG across ~10M variants of 1,400 proteins. "
            "Fast, CPU-friendly, open source (Apache 2.0)."
        ),
        "thermompnn": (
            "**ThermoMPNN** — Graph transformer trained on the Megascale cDNA display "
            "dataset (~770K experimental ΔΔG measurements). State-of-the-art on the "
            "Megascale benchmark. Open source (MIT)."
        ),
        "foldx": (
            "**FoldX** — Empirical force field with terms for van der Waals, electrostatics, "
            "solvation, and backbone strain. Calibrated against directed-evolution data. "
            "Free for academic use; commercial license required otherwise."
        ),
    }
    lines = ["## Methods\n"]
    for method in methods_available:
        if method in descriptions:
            lines.append(f"- {descriptions[method]}")
    return "\n".join(lines)


def _render_references(methods_available: list[str]) -> str:
    lines = ["## References\n"]
    for method in methods_available:
        if method in _CITATIONS:
            lines.append(f"- {_CITATIONS[method]}")
    lines.append(f"- Validation context: {_VALIDATION_CITATION}")
    return "\n".join(lines)


def _render_reproducibility(output_dir: Path) -> str:
    # Files are written by stability_predictor._write_outputs after this report
    # is rendered; we list the artefacts that will be present rather than
    # probing the filesystem mid-render.
    del output_dir  # unused, kept for signature stability
    lines = [
        "## Reproducibility",
        "",
        "This run produces the following replay artefacts:",
        "",
        "- `reproducibility/commands.sh` — exact CLI invocation",
        "- `reproducibility/environment.yml` — method versions and dependencies",
        "- `reproducibility/input_checksum.txt` — SHA-256 of inputs",
        "- `reproducibility/output_checksum.txt` — SHA-256 of all outputs",
        "",
        "**Replay**: `bash reproducibility/commands.sh`",
        "**Verify**: `sha256sum -c reproducibility/output_checksum.txt`",
    ]
    return "\n".join(lines)


def _render_disclaimer() -> str:
    return (
        "## Disclaimer\n\n"
        "This skill provides computational predictions for research planning only. "
        "It does not constitute clinical advice, and predictions should not be used "
        "as the sole basis for clinical, diagnostic, or therapeutic decisions. "
        "Experimental validation is required for any high-stakes use."
    )


# ----------------------------------------------------------------------------
# Small formatters
# ----------------------------------------------------------------------------

def _format_direction(direction: Direction) -> str:
    if direction == Direction.UNKNOWN:
        return "—"
    return f"**{direction.value}**"


def _format_confidence(cr: ConsensusResult) -> str:
    if cr.confidence == Confidence.NONE:
        return "—"
    method_summary = f"{cr.n_methods_succeeded}/{cr.n_methods_attempted} methods"
    return f"{cr.confidence.value} ({method_summary})"
