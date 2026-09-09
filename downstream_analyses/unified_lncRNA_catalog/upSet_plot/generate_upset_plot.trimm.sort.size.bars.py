#!/usr/bin/env python3
"""
generate_upset_plot.py

Stage 2 of the pipeline: takes the consensus_matrix.tsv produced by
build_consensus_upset_matrix.py and renders the final UpSet plot, with the
set-size (bottom-left) bars showing TRUE raw gene counts per catalog rather
than "number of consensus loci this catalog appears in".

Why this exists
----------------
By default, an UpSet plot's set-size bars are just the marginal sum of the
same presence/absence matrix used for the intersection bars. That's correct
for the intersection bars (which should stay at the loci level, so they
don't re-inflate), but it silently undercounts catalogs whose genes were
fragmented and merged into shared loci during clustering (see
fragmentation_index.tsv from stage 1). This script overrides just the
set-size panel with the true counts, which is safe because `upsetplot`
stores that panel as an independent pandas.Series with no back-effect on
the intersection matrix.

IMPORTANT CAVEAT (put this in your figure legend / methods):
Once overridden, a category's total bar no longer equals the sum of the
intersection bars it participates in -- that invariant is intentionally
broken here. A catalog whose genes were heavily fragmented and merged into
shared loci (e.g. NONCODE, MiTranscriptome) will show a total bar taller
than the sum of its intersection bars. State explicitly in the caption that
set-size bars = raw annotated gene counts, while intersection bars =
counts of consensus loci after exon-overlap clustering.

Requirements
------------
pip install pandas upsetplot matplotlib "pandas<2.2" --break-system-packages
(pandas<2.2 works around a copy-on-write bug in upsetplot==0.9.0 as of
writing -- check for a newer upsetplot release that fixes this before
pinning down pandas long-term.)

Usage
-----
    python3 generate_upset_plot.py \
        --matrix results/consensus_matrix.tsv \
        --catalogs GENCODE47 GENCODE27 lncRNAmerge NONCODE MiTranscriptome RefSeq \
        --out results/upset_final.png \
        --title "Consensus lncRNA loci across annotation catalogs"

    # Use loci-present totals instead (the plain/default behavior):
    python3 generate_upset_plot.py --matrix results/consensus_matrix.tsv \
        --catalogs GENCODE47 GENCODE27 --out results/upset_default.png --no-true-totals
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def build_upset(matrix_df: pd.DataFrame, catalogs: list, true_totals: bool,
                 min_subset_size: int, sort_by: str, max_subset_rank: int,
                 max_degree: int, sort_categories_by: str, element_size,
                 totals_plot_elements: int):
    """element_size: None -> fit content into whatever figsize you set on the
    Figure before calling upset.plot(fig=fig) (canvas size is respected
    exactly). A number (upsetplot's own default is 32) -> upsetplot instead
    IGNORES the figure's current size and expands/shrinks the canvas to fit
    that many points per row/column. These two modes are mutually exclusive
    in upsetplot; pass None whenever you need a pinned, exact output size.

    totals_plot_elements: width (in grid units) given to the set-size panel
    on the left. upsetplot's default is 2, which is why those bars come out
    thin and cramped -- 5-8 gives long, prominent bars similar to a
    hand-styled ComplexUpset (R) figure.
    """
    from upsetplot import UpSet, from_indicators

    bool_df = matrix_df[catalogs].astype(bool)
    data = from_indicators(bool_df, data=matrix_df)

    upset = UpSet(
        data,
        subset_size="count",
        show_counts=True,
        sort_by=sort_by,
        sort_categories_by=sort_categories_by,
        min_subset_size=min_subset_size if min_subset_size > 0 else None,
        max_subset_rank=max_subset_rank if max_subset_rank > 0 else None,
        max_degree=max_degree if max_degree > 0 else None,
        facecolor="#4C72B0",
        element_size=element_size,
        totals_plot_elements=totals_plot_elements,
    )

    if true_totals:
        n_cols = [f"{c}_n" for c in catalogs]
        missing = [c for c in n_cols if c not in matrix_df.columns]
        if missing:
            sys.exit(
                f"ERROR: --true-totals requires '*_n' columns from stage 1 "
                f"(build_consensus_upset_matrix.py). Missing: {missing}"
            )
        true_counts = {c: int(matrix_df[f"{c}_n"].sum()) for c in catalogs}
        upset.totals = pd.Series(true_counts).reindex(upset.totals.index)

    return upset


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", type=Path, required=True,
                     help="consensus_matrix.tsv from build_consensus_upset_matrix.py")
    ap.add_argument("--catalogs", nargs="+", required=True,
                     help="Catalog names, matching column names in the matrix (order = plot order)")
    ap.add_argument("--out", type=Path, required=True,
                     help="Output path -- extension decides format (.pdf recommended for many-catalog "
                          "plots: vector, stays legible/editable no matter how many bars; .png also works)")
    ap.add_argument("--title", default="Consensus loci across catalogs")
    ap.add_argument("--no-true-totals", dest="true_totals", action="store_false",
                     help="Disable the true-totals override; use default loci-present totals instead")
    ap.add_argument("--min-subset-size", type=int, default=0,
                     help="Hide intersection bars smaller than this (0 = show all)")
    ap.add_argument("--max-subset-rank", type=int, default=40,
                     help="Show only the top N intersection bars by size (0 = show all -- "
                          "with many catalogs this can mean up to 2^n bars, so avoid 0 for >6 catalogs). "
                          "Default 40 keeps the plot readable.")
    ap.add_argument("--max-degree", type=int, default=0,
                     help="Only show intersections involving at most this many catalogs at once "
                          "(0 = no limit). E.g. --max-degree 4 hides combinations of 5+ catalogs, "
                          "which are usually tiny, numerous, and not the story you're telling.")
    ap.add_argument("--sort-by", default="cardinality", choices=["cardinality", "degree", "-cardinality", "-degree"])
    ap.add_argument("--sort-categories-by", default="input",
                     choices=["input", "-input", "cardinality", "-cardinality"],
                     help="Row order in the matrix panel. 'input' (default) = exactly the order you "
                          "listed in --catalogs, top to bottom. 'cardinality' = largest set first "
                          "(upsetplot's own default, ignores --catalogs order).")
    ap.add_argument("--width", type=float, default=None,
                     help="Exact output canvas width in inches (e.g. 14, to match a reference figure's "
                          "artboard). If set together with --height, overrides the auto-sized figure "
                          "AND disables bbox_inches='tight', so the saved PDF/PNG page is exactly "
                          "width x height -- matching how a fixed Illustrator artboard behaves, rather "
                          "than auto-cropping to whatever the content happens to fill.")
    ap.add_argument("--height", type=float, default=None,
                     help="Exact output canvas height in inches. Must be given together with --width.")
    ap.add_argument("--element-size", type=int, default=None,
                     help="upsetplot's fixed points-per-row/column layout (upsetplot's own default is "
                          "32 if you pass a number here). WARNING: setting this to a number makes "
                          "upsetplot ignore --width/--height and resize the canvas to fit that many "
                          "points per row instead. Leave unset (default) to have content fit exactly "
                          "into --width/--height; only set this if you want auto-expanding canvas "
                          "sizing instead of a pinned output size.")
    ap.add_argument("--dpi", type=int, default=200, help="Only affects raster formats like .png; ignored for .pdf")
    ap.add_argument("--totals-plot-elements", type=int, default=6,
                     help="Width (grid units) given to the set-size panel on the left (upsetplot's own "
                          "default is 2, which makes those bars thin/cramped). Default here is 6 for "
                          "longer, more prominent bars. Increase further for very long catalog names.")
    ap.add_argument("--totals-scale-max", type=float, default=None,
                     help="Fix the set-size axis maximum to this value instead of auto-scaling to this "
                          "run's own largest catalog. Use this to match the scale of a reference figure "
                          "(e.g. --totals-scale-max 250000) so bar lengths are visually comparable "
                          "across figures/panels.")
    ap.add_argument("--totals-bar-height", type=float, default=0.85,
                     help="Thickness of each set-size bar, as a fraction of its row (0-1). upsetplot's "
                          "own default is 0.5 (thin); 0.85 gives thick, expressive bars.")
    ap.add_argument("--totals-axis-style", default="clean", choices=["clean", "default"],
                     help="'clean' (default) strips the small stray axis line/ticks from the set-size "
                          "panel, relying on the per-bar count labels only. 'default' keeps upsetplot's "
                          "own minimal 2-tick axis.")
    args = ap.parse_args()

    if (args.width is None) != (args.height is None):
        sys.exit("ERROR: --width and --height must be given together")
    if args.width is not None and args.element_size is not None:
        print("WARNING: --element-size with a numeric value overrides --width/--height "
              "(upsetplot expands the canvas to fit that many points per row instead of "
              "respecting your pinned size). Drop --element-size to get the exact canvas "
              "you asked for, or drop --width/--height to let it auto-size.", file=sys.stderr)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matrix_df = pd.read_csv(args.matrix, sep="\t")
    missing_cats = [c for c in args.catalogs if c not in matrix_df.columns]
    if missing_cats:
        sys.exit(f"ERROR: these catalogs are not columns in {args.matrix}: {missing_cats}")

    upset = build_upset(matrix_df, args.catalogs, args.true_totals, args.min_subset_size,
                         args.sort_by, args.max_subset_rank, args.max_degree, args.sort_categories_by,
                         args.element_size, args.totals_plot_elements)

    n_bars = len(upset.intersections)

    if args.width is not None:
        fig_w, fig_h = args.width, args.height
        use_tight_bbox = False   # pinned size -> exact canvas, no auto-crop
    else:
        fig_w = max(7, min(0.35 * n_bars + 4, 40))   # cap width so it can't run away
        fig_h = max(5, 0.5 * len(args.catalogs) + 3)
        use_tight_bbox = True    # auto-sized -> crop to content as before

    fig = plt.figure(figsize=(fig_w, fig_h))
    axes = upset.plot(fig=fig)
    plt.suptitle(args.title, y=1.02, fontsize=12)

    totals_ax = axes.get("totals")
    if totals_ax is not None:
        if args.totals_scale_max is not None:
            totals_ax.set_xlim(args.totals_scale_max, 0)

        if args.totals_bar_height != 0.5:
            for patch in totals_ax.patches:
                old_h = patch.get_height()
                new_h = args.totals_bar_height
                patch.set_y(patch.get_y() - (new_h - old_h) / 2)
                patch.set_height(new_h)

        if args.totals_axis_style == "clean":
            totals_ax.set_xticks([])
            for spine in totals_ax.spines.values():
                spine.set_visible(False)
            totals_ax.grid(False)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    if use_tight_bbox:
        fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    else:
        fig.savefig(args.out, dpi=args.dpi)   # exact canvas, no cropping
    print(f"wrote {args.out} ({n_bars} intersection bars shown), "
          f"canvas {fig_w:.2f} x {fig_h:.2f} in ({'exact' if not use_tight_bbox else 'tight-cropped'})",
          file=sys.stderr)
    if args.max_subset_rank > 0:
        print(f"  (capped to top {args.max_subset_rank} intersections by --max-subset-rank; "
              f"use --max-subset-rank 0 to show all)", file=sys.stderr)

    if args.true_totals:
        print(
            "\nNOTE: set-size bars show TRUE raw gene counts (from *_n columns); "
            "intersection bars show consensus-locus counts. A category's total "
            "bar may exceed the sum of its intersection bars where genes were "
            "fragmented and merged during clustering -- state this explicitly "
            "in your figure legend.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
