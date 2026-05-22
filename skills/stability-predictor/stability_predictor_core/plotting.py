"""Plot generation for stability-predictor reports.

Three figures are produced when a run has data to plot:

  figures/ddg_heatmap.png        — mutation (rows) × method (cols), colour = ΔΔG
  figures/method_agreement.png   — pairwise scatter of method-vs-method ΔΔG
  figures/per_mutation_bars.png  — grouped bars, one cluster per mutation

matplotlib is the only external dependency. If it is unavailable, plotting
is skipped silently — the report renderer already conditionally includes
figure references based on file presence.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .consensus import ConsensusResult


def render_figures(
    consensus_results: list["ConsensusResult"],
    methods_available: list[str],
    output_dir: Path,
) -> list[Path]:
    """Render the figure set for the report. Returns the list of files written.

    Skips any figure whose data is degenerate (e.g., a scatter with <2
    overlapping predictions). Never raises; logs and continues on failure.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.info("matplotlib not installed; skipping figure generation")
        return []

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    table = _build_table(consensus_results, methods_available)
    if not table["mutations"]:
        return written

    for fn in (_plot_heatmap, _plot_agreement, _plot_per_mutation_bars):
        try:
            path = fn(plt, table, figures_dir)
            if path is not None:
                written.append(path)
        except Exception as exc:  # noqa: BLE001 - plotting is best-effort
            logger.warning("Figure %s failed: %s", fn.__name__, exc)

    return written


def _build_table(
    consensus_results: list["ConsensusResult"],
    methods_available: list[str],
) -> dict:
    """Pivot consensus results into a {mutation: {method: ddg}} table."""
    mutations: list[str] = []
    rows: list[dict[str, float | None]] = []
    for cr in consensus_results:
        mutations.append(str(cr.mutation))
        by_method: dict[str, float | None] = {m: None for m in methods_available}
        for pred in cr.per_method:
            if pred.method in by_method and pred.succeeded:
                by_method[pred.method] = pred.ddg
        rows.append(by_method)
    return {"mutations": mutations, "methods": methods_available, "rows": rows}


# ----------------------------------------------------------------------------
# Individual plots
# ----------------------------------------------------------------------------

def _plot_heatmap(plt, table: dict, figures_dir: Path) -> Path | None:
    """Heatmap of ΔΔG with mutations on rows, methods on columns."""
    mutations = table["mutations"]
    methods = table["methods"]
    if not methods:
        return None

    # Build a 2D array with NaN for missing values
    import numpy as np

    matrix = np.full((len(mutations), len(methods)), np.nan)
    for i, row in enumerate(table["rows"]):
        for j, method in enumerate(methods):
            value = row.get(method)
            if value is not None:
                matrix[i, j] = value

    fig_height = max(2.5, 0.4 * len(mutations) + 1.5)
    fig, ax = plt.subplots(figsize=(max(3.5, 1.2 * len(methods) + 1.5), fig_height))

    vmax = max(abs(np.nanmin(matrix)) if np.any(~np.isnan(matrix)) else 1.0,
               abs(np.nanmax(matrix)) if np.any(~np.isnan(matrix)) else 1.0,
               1.0)
    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods, rotation=30, ha="right")
    ax.set_yticks(range(len(mutations)))
    ax.set_yticklabels(mutations)
    ax.set_title("ΔΔG (kcal/mol) by method")

    # Annotate cells with the numeric value (or '—' for missing)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if np.isnan(value):
                text = "—"
                colour = "#888"
            else:
                text = f"{value:+.2f}"
                colour = "black" if abs(value) < 0.6 * vmax else "white"
            ax.text(j, i, text, ha="center", va="center", color=colour, fontsize=8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("ΔΔG (kcal/mol)")
    fig.tight_layout()

    out = figures_dir / "ddg_heatmap.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def _plot_agreement(plt, table: dict, figures_dir: Path) -> Path | None:
    """Scatter of method-vs-method ΔΔG for the first two methods with data.

    For runs with only one method, this figure is skipped.
    """
    methods = table["methods"]
    if len(methods) < 2:
        return None

    # Find the first two methods that have at least 2 overlapping predictions
    import numpy as np

    pair = None
    for i in range(len(methods)):
        for j in range(i + 1, len(methods)):
            xs, ys = [], []
            for row in table["rows"]:
                xi = row.get(methods[i])
                yi = row.get(methods[j])
                if xi is not None and yi is not None:
                    xs.append(xi)
                    ys.append(yi)
            if len(xs) >= 2:
                pair = (methods[i], methods[j], xs, ys)
                break
        if pair:
            break

    if pair is None:
        return None

    m_x, m_y, xs, ys = pair
    xs_arr = np.array(xs)
    ys_arr = np.array(ys)

    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.scatter(xs_arr, ys_arr, s=42, edgecolors="black", linewidths=0.6, alpha=0.85)

    lim = max(float(np.max(np.abs(xs_arr))), float(np.max(np.abs(ys_arr))), 1.0) * 1.15
    ax.plot([-lim, lim], [-lim, lim], linestyle="--", color="grey", linewidth=1.0,
            label="y = x")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(f"{m_x} ΔΔG (kcal/mol)")
    ax.set_ylabel(f"{m_y} ΔΔG (kcal/mol)")
    ax.set_title(f"Method agreement: {m_x} vs {m_y}")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    out = figures_dir / "method_agreement.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def _plot_per_mutation_bars(plt, table: dict, figures_dir: Path) -> Path | None:
    """Grouped bars per mutation, one bar per method."""
    mutations = table["mutations"]
    methods = table["methods"]
    if not mutations or not methods:
        return None

    import numpy as np

    n_mut = len(mutations)
    n_meth = len(methods)
    width = 0.8 / n_meth
    indices = np.arange(n_mut)

    fig_width = max(5.0, 0.7 * n_mut + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, 3.5))

    for k, method in enumerate(methods):
        values = []
        for row in table["rows"]:
            v = row.get(method)
            values.append(v if v is not None else 0.0)
        offsets = indices + (k - (n_meth - 1) / 2) * width
        ax.bar(offsets, values, width=width, label=method,
               edgecolor="black", linewidth=0.4)

    ax.axhline(0, color="black", linewidth=0.6)
    ax.axhline(1.0, color="red", linewidth=0.5, linestyle=":")
    ax.axhline(-1.0, color="blue", linewidth=0.5, linestyle=":")
    ax.set_xticks(indices)
    ax.set_xticklabels(mutations, rotation=30, ha="right")
    ax.set_ylabel("ΔΔG (kcal/mol)")
    ax.set_title("Per-mutation ΔΔG by method")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()

    out = figures_dir / "per_mutation_bars.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out