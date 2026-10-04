import os, glob
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm, ListedColormap
from matplotlib.patches import Patch

from make_heat_map import (
    MUTATED_FILES_DIR,
    MUTATED_PATTERN,
    MAX_PATTERN_LEN,
    HEATMAP_DIR,
    NODATA_COLOR,
    ZERO_COLOR,
    count_mutations_per_cell,
    read_callable_per_cell,
    sample_name,
    count_file_for,
)

# Sample name prefix (case-insensitive) identifying the HG005 sample.
HG005_SAMPLE_PREFIX = "hg005"


def _labels(per_locus: bool):
    """(colorbar label, title metric) for the chosen denominator."""
    if per_locus:
        return "mutations per callable locus", "mutation density (mutations / callable locus)"
    return "mutations per callable base", "mutation density (mutations / callable base)"


def _mode_suffix(per_locus: bool) -> str:
    return "_per_locus" if per_locus else "_per_base"


def _compute_combined_density(exclude_hg005: bool = False, per_locus: bool = False):
    """Pool mutations and the callable denominator across all GIAB samples and
    return (density, max_repeats, samples_used).

        density(pattern_length, num_repeats) = total_mutations / total_denominator

    The denominator is summed over all included samples per (pattern_length,
    num_repeats) cell and is either:
      * per_locus=False (default): total callable *bases*
                                   (pattern_length * num_repeats * num_callable_loci)
      * per_locus=True:            total callable *loci* (num_callable_loci)

    When `exclude_hg005` is True the HG005 sample is left out. NaN cells mean the
    denominator is zero (density undefined).
    """
    total_mutations = defaultdict(int)   # summed mutations per (pattern_length, num_repeats)
    total_denominator = defaultdict(int)  # summed callable loci or bases per cell
    samples_used = []

    for mutated_fp in sorted(glob.glob(os.path.join(MUTATED_FILES_DIR, MUTATED_PATTERN))):
        name = sample_name(mutated_fp)
        if exclude_hg005 and name.lower().startswith(HG005_SAMPLE_PREFIX):
            continue
        count_fp = count_file_for(mutated_fp)
        if not count_fp:
            print(f"{name}: no matching count file, skipping")
            continue

        samples_used.append(name)
        mutations = count_mutations_per_cell(mutated_fp)
        callable_loci = read_callable_per_cell(count_fp)

        for (pattern_len, num_repeats), num_callable in callable_loci.items():
            if not (1 <= pattern_len <= MAX_PATTERN_LEN) or num_repeats < 1 or num_callable <= 0:
                continue
            if per_locus:
                total_denominator[(pattern_len, num_repeats)] += num_callable
            else:
                total_denominator[(pattern_len, num_repeats)] += pattern_len * num_repeats * num_callable
        for key, muts in mutations.items():
            total_mutations[key] += muts

    max_repeats = max((nr for _, nr in total_denominator), default=0)
    density = np.full((MAX_PATTERN_LEN, max_repeats), np.nan)
    placed = 0
    for (pattern_len, num_repeats), denom in total_denominator.items():
        if denom <= 0:
            continue
        muts = total_mutations.get((pattern_len, num_repeats), 0)
        density[pattern_len - 1, num_repeats - 1] = muts / denom
        placed += muts

    total_muts = sum(total_mutations.values())
    if placed < total_muts:
        print(f"  note: {total_muts - placed}/{total_muts} mutations fell on (pattern_length, num_repeats) "
              f"cells with no callable count and were not placed")

    return density, max_repeats, samples_used


def _pad_cols(density: np.ndarray, target_cols: int) -> np.ndarray:
    """Right-pad a density matrix with NaN columns so it has `target_cols` columns
    (aligns the num_repeats axis of two matrices for side-by-side comparison)."""
    rows, cols = density.shape
    if cols >= target_cols:
        return density
    return np.hstack([density, np.full((rows, target_cols - cols), np.nan)])


def _log_norm(*matrices):
    """Shared LogNorm spanning the positive densities of the given matrices, or None."""
    positive = np.concatenate([m[m > 0].ravel() for m in matrices]) if matrices else np.array([])
    if not positive.size:
        return None
    vmin, vmax = positive.min(), positive.max()
    if vmin == vmax:
        vmin = vmax / 10.0
    return LogNorm(vmin=vmin, vmax=vmax)


def _draw_panel(ax, density: np.ndarray, max_repeats: int, norm):
    """Draw one heat map on `ax` using a caller-supplied (shared) color norm, and
    return the positive-density image handle (for a shared colorbar)."""
    # Layer 1 (background): callable-but-zero-mutation cells vs no-callable cells.
    callable_bg = np.ma.masked_invalid(np.where(np.isfinite(density), 0.0, np.nan))
    zero_cmap = ListedColormap([ZERO_COLOR])
    zero_cmap.set_bad(NODATA_COLOR)
    ax.imshow(callable_bg, aspect="auto", origin="lower", cmap=zero_cmap, vmin=0, vmax=1)

    # Layer 2 (overlay): positive densities on the shared log scale.
    positive = np.ma.masked_where(~(density > 0), density)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(alpha=0.0)
    im = ax.imshow(positive, aspect="auto", origin="lower", cmap=cmap, norm=norm)

    ax.set_xlabel("number of repeats")
    ax.set_ylabel("pattern length")
    ax.set_yticks(range(MAX_PATTERN_LEN))
    ax.set_yticklabels(range(1, MAX_PATTERN_LEN + 1))
    step = max(1, max_repeats // 20)
    ax.set_xticks(range(0, max_repeats, step))
    ax.set_xticklabels(range(1, max_repeats + 1, step))
    return im


def _legend_handles():
    return [Patch(facecolor=ZERO_COLOR, edgecolor="gray", label="callable, 0 mutations"),
            Patch(facecolor=NODATA_COLOR, edgecolor="gray", label="no callable loci")]


def combined_heat_map(exclude_hg005: bool = False, per_locus: bool = False) -> np.ndarray:
    """Combined mutation-density heat map across all GIAB samples (one panel).

    When `exclude_hg005` is True the HG005 sample is left out (produces a heat map
    for all the other samples); otherwise every sample is included. `per_locus`
    selects the denominator: callable loci (True) or callable bases (False). Saves a
    PNG and returns the density matrix (rows: pattern_length 1..15, cols: num_repeats).
    """
    density, max_repeats, samples_used = _compute_combined_density(exclude_hg005, per_locus)
    name = f"combined_{'excluding_HG005' if exclude_hg005 else 'all_samples'}{_mode_suffix(per_locus)}"
    print(f"{name}: combined heat map over {len(samples_used)} samples")

    os.makedirs(HEATMAP_DIR, exist_ok=True)
    cbar_label, title_metric = _labels(per_locus)
    norm = _log_norm(density)
    fig, ax = plt.subplots(figsize=(max(6, max_repeats * 0.3), 6), constrained_layout=True)
    im = _draw_panel(ax, density, max_repeats, norm)
    ax.set_title(f"{name}\n{title_metric}")
    if norm is not None:
        fig.colorbar(im, ax=ax, label=cbar_label)
    ax.legend(handles=_legend_handles(), loc="upper right", framealpha=0.9, fontsize=8)

    out_fp = os.path.join(HEATMAP_DIR, f"{name}.heatmap.png")
    fig.savefig(out_fp, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_fp}")
    return density


def combined_heat_maps_side_by_side(per_locus: bool = False) -> str:
    """Render the two combined heat maps (all samples vs excluding HG005) as two
    panels on a single figure, sharing one log color scale and a common
    num_repeats axis so they can be compared directly. `per_locus` selects the
    denominator (callable loci vs callable bases). Saves a PNG and returns its path.
    """
    density_all, max_all, samples_all = _compute_combined_density(exclude_hg005=False, per_locus=per_locus)
    density_excl, max_excl, samples_excl = _compute_combined_density(exclude_hg005=True, per_locus=per_locus)

    common_max = max(max_all, max_excl)
    density_all = _pad_cols(density_all, common_max)
    density_excl = _pad_cols(density_excl, common_max)

    # Shared log color scale spanning positive densities of BOTH panels.
    norm = _log_norm(density_all, density_excl)
    cbar_label, title_metric = _labels(per_locus)

    os.makedirs(HEATMAP_DIR, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(max(12, common_max * 0.6), 6),
                             sharey=True, constrained_layout=True)
    im = _draw_panel(axes[0], density_all, common_max, norm)
    _draw_panel(axes[1], density_excl, common_max, norm)
    axes[0].set_title(f"all samples (n={len(samples_all)})")
    axes[1].set_title(f"excluding HG005 (n={len(samples_excl)})")
    fig.suptitle(title_metric)

    if norm is not None:
        fig.colorbar(im, ax=axes, label=cbar_label)
    axes[1].legend(handles=_legend_handles(), loc="upper right", framealpha=0.9, fontsize=8)

    out_fp = os.path.join(HEATMAP_DIR, f"combined_side_by_side{_mode_suffix(per_locus)}.heatmap.png")
    fig.savefig(out_fp, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_fp}")
    return out_fp


# ---------------------------------------------------------------------------
# Publication figure
# ---------------------------------------------------------------------------
# Standalone, submission-ready rendering of the pooled density matrix: no title
# (the journal caption carries the description), no sample names, real vector
# text in the PDF, and type sizes that stay legible at the printed width.

FIG_WIDTH_IN = 180 / 25.4  # 180 mm, the usual double-column width

# Lighter than the exploratory NODATA_COLOR / ZERO_COLOR: the measured-zero
# region covers most of the panel, and at print size a mid-grey there competes
# with the dark end of viridis. Substitute the shared constants to match the
# other figures in this module.
PUB_ZERO_COLOR = "#dcdcdc"    # callable, zero mutations observed
PUB_NODATA_COLOR = "#ffffff"  # no callable loci of this (pattern_length, num_repeats)

PUB_RCPARAMS = {
    "font.family": "sans-serif",
    # Liberation Sans / Nimbus Sans are metric-compatible with Arial / Helvetica
    # and ship on most Linux boxes; DejaVu (matplotlib's default) is neither, and
    # is noticeably wider. First installed name wins.
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
    "pdf.fonttype": 42,   # embed TrueType so text stays editable/selectable
    "ps.fonttype": 42,
    "font.size": 7,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "axes.labelpad": 3.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "legend.frameon": False,
    "legend.handlelength": 1.1,
    "legend.handleheight": 1.1,
    "legend.columnspacing": 1.6,
}


# Margins in inches, laid out explicitly rather than by constrained_layout so the
# printed panel size is reproducible and the colorbar height tracks the panel
# exactly. Tuned for the type sizes in PUB_RCPARAMS.
PUB_MARGIN_LEFT = 0.62     # y-axis label + two-digit tick labels
PUB_MARGIN_RIGHT = 0.64    # colorbar gap + bar + colorbar label
PUB_MARGIN_BOTTOM = 0.64   # x tick labels + x-axis label + legend row
PUB_MARGIN_TOP = 0.06
PUB_CBAR_GAP = 0.08
PUB_CBAR_WIDTH = 0.10


def _figure_layout(max_repeats: int):
    """(figsize, panel rect, colorbar rect) for a 180 mm-wide panel of
    MAX_PATTERN_LEN rows by `max_repeats` columns. Rects are in figure fractions.

    Cells are deliberately taller than they are wide: with ~75 repeat columns a
    square cell leaves the pattern-length tick labels overlapping.
    """
    axes_w = FIG_WIDTH_IN - PUB_MARGIN_LEFT - PUB_MARGIN_RIGHT
    col_w = axes_w / max(max_repeats, 1)
    row_h = min(0.20, max(0.11, col_w * 1.6))
    axes_h = row_h * MAX_PATTERN_LEN

    width = FIG_WIDTH_IN
    height = axes_h + PUB_MARGIN_BOTTOM + PUB_MARGIN_TOP
    panel = [PUB_MARGIN_LEFT / width, PUB_MARGIN_BOTTOM / height,
             axes_w / width, axes_h / height]
    cbar = [(PUB_MARGIN_LEFT + axes_w + PUB_CBAR_GAP) / width, PUB_MARGIN_BOTTOM / height,
            PUB_CBAR_WIDTH / width, axes_h / height]
    return (width, height), panel, cbar


def _tick_step(n: int, max_ticks: int = 16) -> int:
    """Smallest 1/2/5/10/20/25/50-style step keeping the tick count under `max_ticks`."""
    for step in (1, 2, 5, 10, 20, 25, 50, 100):
        if n / step <= max_ticks:
            return step
    return max(1, n // max_ticks)


def _draw_publication_panel(ax, density: np.ndarray, max_repeats: int, norm):
    """Same two-layer construction as `_draw_panel`, restyled for print."""
    # Layer 1 (background): measured zeros vs cells with no callable loci.
    callable_bg = np.ma.masked_invalid(np.where(np.isfinite(density), 0.0, np.nan))
    zero_cmap = ListedColormap([PUB_ZERO_COLOR])
    zero_cmap.set_bad(PUB_NODATA_COLOR)
    ax.imshow(callable_bg, aspect="auto", origin="lower", cmap=zero_cmap,
              vmin=0, vmax=1, interpolation="nearest")

    # Layer 2 (overlay): positive densities on the log scale.
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(alpha=0.0)  # let the background layer show through
    im = ax.imshow(np.ma.masked_where(~(density > 0), density), aspect="auto",
                   origin="lower", cmap=cmap, norm=norm, interpolation="nearest")

    ax.set_xlabel("Number of repeats")
    ax.set_ylabel("Pattern length (bp)")

    ystep = _tick_step(MAX_PATTERN_LEN)
    yticks = list(range(0, MAX_PATTERN_LEN, ystep))
    ax.set_yticks(yticks)
    ax.set_yticklabels([t + 1 for t in yticks])

    xstep = _tick_step(max_repeats)
    # Always anchor the first tick at 1 repeat, then step from there.
    xticks = [0] + [r - 1 for r in range(xstep, max_repeats + 1, xstep)]
    ax.set_xticks(xticks)
    ax.set_xticklabels([t + 1 for t in xticks])

    ax.tick_params(top=False, right=False)
    return im


def publication_heat_map(per_locus: bool = True, exclude_hg005: bool = True,
                         basename: str = "") -> str:
    """Submission-ready combined mutation-density heat map.

    Pools mutations and the callable denominator across the GIAB samples and
    renders a single unlabelled panel as vector PDF plus a 600 dpi PNG. Neither
    the file name nor the figure mentions which samples were pooled. Returns the
    PDF path.
    """
    density, max_repeats, samples_used = _compute_combined_density(exclude_hg005, per_locus)
    if max_repeats == 0:
        raise RuntimeError("no callable loci found; nothing to plot")

    basename = basename or f"mutation_density{_mode_suffix(per_locus)}"
    cbar_label, _ = _labels(per_locus)
    cbar_label = cbar_label[0].upper() + cbar_label[1:]
    norm = _log_norm(density)
    figsize, panel_rect, cbar_rect = _figure_layout(max_repeats)

    os.makedirs(HEATMAP_DIR, exist_ok=True)
    with plt.rc_context(PUB_RCPARAMS):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes(panel_rect)
        im = _draw_publication_panel(ax, density, max_repeats, norm)

        if norm is not None:
            cbar = fig.colorbar(im, cax=fig.add_axes(cbar_rect))
            cbar.set_label(cbar_label)
            cbar.ax.tick_params(labelsize=7, width=0.6, length=2.5)
            cbar.ax.tick_params(which="minor", width=0.4, length=1.5)
            cbar.ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation(base=10))
            cbar.ax.yaxis.set_minor_formatter(ticker.NullFormatter())
            cbar.outline.set_linewidth(0.6)

        handles = [
            Patch(facecolor=PUB_ZERO_COLOR, edgecolor="#4d4d4d", linewidth=0.5,
                  label="Callable, no mutations observed"),
            Patch(facecolor=PUB_NODATA_COLOR, edgecolor="#4d4d4d", linewidth=0.5,
                  label="No callable loci"),
        ]
        # Centred under the panel (not the whole figure, which includes the colorbar).
        fig.legend(handles=handles, ncol=2, loc="lower center",
                   bbox_to_anchor=(panel_rect[0] + panel_rect[2] / 2, 0.005),
                   bbox_transform=fig.transFigure)

        pdf_fp = os.path.join(HEATMAP_DIR, f"{basename}.pdf")
        png_fp = os.path.join(HEATMAP_DIR, f"{basename}.png")
        fig.savefig(pdf_fp)              # vector, text preserved
        fig.savefig(png_fp, dpi=600)     # raster preview at submission resolution
        plt.close(fig)

    print(f"{basename}: publication figure over {len(samples_used)} samples")
    print(f"  wrote {pdf_fp}")
    print(f"  wrote {png_fp}")
    return pdf_fp


def main():
    # per callable base (original) and per callable locus (new)
    # for per_locus in (False, True):
    #     combined_heat_map(exclude_hg005=False, per_locus=per_locus)
    #     combined_heat_map(exclude_hg005=True, per_locus=per_locus)
    #     combined_heat_maps_side_by_side(per_locus=per_locus)
    publication_heat_map(per_locus=True)


if __name__ == "__main__":
    main()
