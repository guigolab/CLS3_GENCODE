#!/usr/bin/env python3
"""
build_consensus_upset_matrix.py

Builds a gene-level consensus set across multiple lncRNA catalogs using an
"overlap-fraction on exon footprints" rule, to feed an UpSet plot with
non-inflated set sizes.

Pipeline
--------
1. For each catalog, collapse each gene's exons (across all its isoforms)
   into a non-redundant "exon footprint" (merged intervals per gene_id).
2. For every pair of catalogs, run `bedtools intersect -wo -s` between their
   exon footprints and sum shared bp per (geneA, geneB) pair.
3. Convert summed overlap bp into a fraction of the SMALLER gene's total
   footprint length; keep pairs clearing a threshold as graph edges.
4. Take connected components of the resulting graph as unified "consensus
   loci" (this is the multi-catalog analogue of what you already did for
   lncRNA-merge, but done for ALL catalogs simultaneously, not pairwise
   against GENCODE only).
5. Output a presence/absence matrix (consensus_locus x catalog) ready for
   UpSetR / ComplexUpset (R) or the `upsetplot` package (Python).

Requirements
------------
- bedtools on PATH
- pandas, networkx  (pip install pandas networkx)

Input format
------------
One BED6 file per catalog, exon-level, already tagged with gene_id in
column 4, e.g.:

    chr1  11869  12227  ENSG00000223972  .  +
    chr1  12613  12721  ENSG00000223972  .  +
    ...

Usage
-----
    python3 build_consensus_upset_matrix.py \
        --catalogs GENCODE47=gencode47.exons.bed GENCODE27=gencode27.exons.bed \
                   lncRNAmerge=merge.exons.bed NONCODE=noncode.exons.bed \
        --min-frac 0.3 \
        --outdir results/

Notes
-----
- Strand-aware throughout (-s). Change to unstranded by editing INTERSECT_ARGS
  if needed for any catalog.
- min-frac is applied to the SMALLER gene's footprint (not reciprocal), which
  is deliberately tolerant of fragment-vs-full-length mismatches. Switch to
  reciprocal by editing `frac_smaller` -> `frac_reciprocal` in build_edges().
- This treats footprint overlap for single-exon genes with the same rule.
  If you want a stricter rule specifically for single-exon genes, filter
  them out beforehand and handle separately (see discussion in the paper's
  methods notes).
"""

import argparse
import itertools
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import networkx as nx

BED_COLS = ["chrom", "start", "end", "gene_id", "score", "strand"]


# ---------------------------------------------------------------------------
# Step 1: per-gene exon footprint (merge exons within a gene, not across genes)
# ---------------------------------------------------------------------------
def merge_intervals(intervals):
    """intervals: list of (start, end), sorted or not. Returns merged, sorted list."""
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def build_exon_footprints(bed_path: Path) -> pd.DataFrame:
    """
    Read an exon-level BED6, merge exons per (gene_id, chrom, strand),
    return a long-format DataFrame of footprint blocks:
    chrom, start, end, gene_id, score, strand
    """
    df = pd.read_csv(bed_path, sep="\t", header=None, names=BED_COLS)

    rows = []
    for (gene_id, chrom, strand), sub in df.groupby(["gene_id", "chrom", "strand"]):
        merged = merge_intervals(list(zip(sub["start"], sub["end"])))
        for s, e in merged:
            rows.append((chrom, s, e, gene_id, ".", strand))

    out = pd.DataFrame(rows, columns=BED_COLS)
    return out


def footprint_lengths(footprint_df: pd.DataFrame) -> pd.Series:
    """Total exonic bp per gene_id (sum of merged block lengths)."""
    lengths = (footprint_df["end"] - footprint_df["start"]).groupby(
        footprint_df["gene_id"]
    ).sum()
    return lengths


def write_bed(df: pd.DataFrame, path: Path):
    df.sort_values(["chrom", "start", "end"]).to_csv(
        path, sep="\t", header=False, index=False
    )


# ---------------------------------------------------------------------------
# Step 2: pairwise overlap via bedtools intersect -wo
# ---------------------------------------------------------------------------
INTERSECT_ARGS = ["-s", "-wo"]  # strand-aware, report overlap bp as last column


def pairwise_overlap_bp(bed_a: Path, bed_b: Path) -> pd.DataFrame:
    """
    Run bedtools intersect -wo between two footprint BEDs, aggregate summed
    overlap bp per (gene_a, gene_b) pair (sums across all their exon blocks).
    """
    cmd = ["bedtools", "intersect", "-a", str(bed_a), "-b", str(bed_b)] + INTERSECT_ARGS
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    if not result.stdout.strip():
        return pd.DataFrame(columns=["gene_a", "gene_b", "overlap_bp"])

    # -a columns (6) + -b columns (6) + overlap_bp (1) = 13 columns
    cols = [f"a_{c}" for c in BED_COLS] + [f"b_{c}" for c in BED_COLS] + ["overlap_bp"]
    df = pd.read_csv(
        __import__("io").StringIO(result.stdout), sep="\t", header=None, names=cols
    )

    agg = (
        df.groupby(["a_gene_id", "b_gene_id"])["overlap_bp"]
        .sum()
        .reset_index()
        .rename(columns={"a_gene_id": "gene_a", "b_gene_id": "gene_b"})
    )
    return agg


# ---------------------------------------------------------------------------
# Step 3: convert to overlap fraction, threshold, build edge list
# ---------------------------------------------------------------------------
def build_edges(overlap_df, lengths_a, lengths_b, min_frac):
    """
    fraction = overlap_bp / min(len(gene_a), len(gene_b))
    i.e. fraction of the SMALLER gene's footprint that is covered.
    """
    df = overlap_df.copy()
    df["len_a"] = df["gene_a"].map(lengths_a)
    df["len_b"] = df["gene_b"].map(lengths_b)
    df["frac_smaller"] = df["overlap_bp"] / df[["len_a", "len_b"]].min(axis=1)
    # (optional stricter alternative, reciprocal on both):
    # df["frac_reciprocal"] = df[["overlap_bp"]].div(df[["len_a","len_b"]].max(axis=1))
    return df[df["frac_smaller"] >= min_frac][["gene_a", "gene_b", "frac_smaller"]]


# ---------------------------------------------------------------------------
# Step 4 + 5: connected components -> presence/absence matrix
# ---------------------------------------------------------------------------
def build_consensus_matrix(catalog_gene_ids: dict, edges_by_pair: dict):
    """
    catalog_gene_ids: {catalog_name: set(gene_ids)}
    edges_by_pair:    {(cat_a, cat_b): edge_df with columns gene_a, gene_b}

    Returns
    -------
    matrix_df : DataFrame for the UpSet plot itself.
        One row per consensus locus. For each catalog there are TWO columns:
          - "<cat>"    : 0/1 presence (this is what UpSetR/upsetplot consumes)
          - "<cat>_n"  : how many raw genes from that catalog fell into this
                         locus. This is the column that reveals 1-to-many
                         splits (e.g. one full-length gene in catalog A
                         matching two fragmentary/split genes in catalog B) --
                         a >1 value here means genes are being merged into
                         the locus, not lost. Presence alone can't show this,
                         since 2 genes and 1 gene from the same catalog both
                         just read as "1" in the presence column.
    components  : list of sets of node ids (kept for diagnostics)
    node_catalog: dict node_id -> catalog name (kept for diagnostics)
    """
    G = nx.Graph()

    node_catalog = {}
    node_gene = {}
    for cat, genes in catalog_gene_ids.items():
        for g in genes:
            node_id = f"{cat}::{g}"  # namespace gene ids by catalog to avoid collisions
            G.add_node(node_id)
            node_catalog[node_id] = cat
            node_gene[node_id] = g

    for (cat_a, cat_b), edge_df in edges_by_pair.items():
        for _, row in edge_df.iterrows():
            G.add_edge(f"{cat_a}::{row['gene_a']}", f"{cat_b}::{row['gene_b']}")

    components = list(nx.connected_components(G))

    matrix_rows = []
    for i, comp in enumerate(components):
        locus_id = f"locus_{i:06d}"
        counts = {cat: 0 for cat in catalog_gene_ids}
        for node in comp:
            counts[node_catalog[node]] += 1

        row = {"consensus_locus": locus_id, "n_members": len(comp)}
        for cat, n in counts.items():
            row[cat] = 1 if n > 0 else 0   # presence column, for UpSet
            row[f"{cat}_n"] = n            # raw gene count column, for auditing
        matrix_rows.append(row)

    matrix_df = pd.DataFrame(matrix_rows)
    return matrix_df, components, node_catalog, node_gene


def detect_split_events(matrix_df: pd.DataFrame, catalog_names: list) -> pd.DataFrame:
    """
    Loci where >1 gene from the SAME catalog collapsed into one consensus
    locus (the "1 gene in X matches 2+ genes in Y" case). Not necessarily an
    error -- can be a genuine biological split/merge between annotations --
    but always worth a look, and this is the table that shows you the raw
    gene identities that presence/absence alone would hide.
    """
    n_cols = [f"{c}_n" for c in catalog_names]
    has_split = (matrix_df[n_cols] > 1).any(axis=1)
    cols = ["consensus_locus", "n_members"] + \
        [c for pair in zip(catalog_names, n_cols) for c in pair]
    return matrix_df.loc[has_split, cols].sort_values("n_members", ascending=False)


def dump_components_long(components, node_catalog, node_gene) -> pd.DataFrame:
    """Long-format table: consensus_locus, catalog, gene_id -- for full audit."""
    rows = []
    for i, comp in enumerate(components):
        locus_id = f"locus_{i:06d}"
        for node in comp:
            rows.append((locus_id, node_catalog[node], node_gene[node]))
    return pd.DataFrame(rows, columns=["consensus_locus", "catalog", "gene_id"])


def add_relationship_type(matrix_df: pd.DataFrame, catalog_names: list) -> pd.DataFrame:
    """
    Label each locus by its cardinality pattern across catalogs, using the
    *_n (raw gene count) columns already in matrix_df:
      - "empty"     : 0 catalogs present (shouldn't occur, sanity guard)
      - "singleton" : only 1 catalog present at all (no comparison possible)
      - "1:1"       : every present catalog contributes exactly 1 gene
      - "1:many"    : present catalogs contribute a mix of 1s and >1s
      - "many:many" : every present catalog contributes >1 gene
    This is what lets you report split events as a proportion rather than
    trying to force them into the UpSet bar heights themselves.
    """
    df = matrix_df.copy()
    n_cols = [f"{c}_n" for c in catalog_names]

    def classify(row):
        present_counts = [row[c] for c in n_cols if row[c] > 0]
        if len(present_counts) == 0:
            return "empty"
        if len(present_counts) == 1:
            return "singleton"
        if all(c == 1 for c in present_counts):
            return "1:1"
        if all(c > 1 for c in present_counts):
            return "many:many"
        return "1:many"

    df["relationship_type"] = df.apply(classify, axis=1)
    return df


def fragmentation_index(matrix_df: pd.DataFrame, catalog_gene_ids: dict) -> pd.DataFrame:
    """
    Per catalog: raw_gene_count vs n_loci_present (the number this catalog
    contributes to on the UpSet set-size side bar), and their ratio.

    fragmentation_index = raw_gene_count / n_loci_present

    ~1.0  -> this catalog's genes show up ~1-per-locus; the UpSet bar for
             this catalog faithfully represents its true gene count.
    >>1.0 -> this catalog is having multiple genes absorbed into shared
             loci; its UpSet set-size bar is an UNDERCOUNT of its true
             gene contribution by roughly that factor. Catalogs with
             fragmentary transcript models (MiTranscriptome, NONCODE, etc.)
             are expected to show the highest values here.
    """
    rows = []
    for cat, genes in catalog_gene_ids.items():
        raw = len(genes)
        n_loci = int(matrix_df[cat].sum())
        n_true = int(matrix_df[f"{cat}_n"].sum())
        frag_idx = raw / n_loci if n_loci > 0 else float("nan")
        rows.append({
            "catalog": cat,
            "raw_gene_count": raw,
            "n_loci_present": n_loci,
            "genes_recovered_check": n_true,  # should equal raw_gene_count
            "fragmentation_index": round(frag_idx, 3),
        })
    return pd.DataFrame(rows).sort_values("fragmentation_index", ascending=False)


def plot_dual_bar(frag_df: pd.DataFrame, out_path: Path):
    """
    Companion bar chart: raw gene count vs. loci-present count, per catalog,
    side by side. This is the visual answer to "did the UpSet plot catch all
    of this catalog's genes" -- the gap between the two bars for a given
    catalog IS the fragmentation effect; the UpSet plot itself cannot show
    this (see fragmentation_index docstring for why).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    df = frag_df.sort_values("fragmentation_index", ascending=False)
    x = np.arange(len(df))
    width = 0.38

    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(df)), 5))
    ax.bar(x - width / 2, df["raw_gene_count"], width, label="Raw gene count", color="#1b9e77")
    ax.bar(x + width / 2, df["n_loci_present"], width, label="Loci present (UpSet set size)", color="#d95f02")

    ax.set_xticks(x)
    ax.set_xticklabels(df["catalog"], rotation=45, ha="right")
    ax.set_ylabel("Count")
    ax.set_title("Raw gene count vs. UpSet set size, per catalog")
    ax.legend()

    for i, row in enumerate(df.itertuples()):
        ax.annotate(f"{row.fragmentation_index}x", (i, max(row.raw_gene_count, row.n_loci_present)),
                    textcoords="offset points", xytext=(0, 4), ha="center", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Diagnostics: flag likely chaining artifacts
# ---------------------------------------------------------------------------
def flag_outlier_components(matrix_df: pd.DataFrame, catalog_names: list, pct=99):
    """
    Components with unusually many member genes or unusually many catalogs
    are worth manual review -- likely bridge/chaining artifacts even after
    switching to exon-footprint overlap. (n_members counts every raw gene
    absorbed, so a component built from a 1-to-many split already shows up
    here as larger than a simple 1-to-1 match.)
    """
    thresh = matrix_df["n_members"].quantile(pct / 100)
    n_catalogs = matrix_df[catalog_names].sum(axis=1)
    flagged = matrix_df[(matrix_df["n_members"] >= thresh) | (n_catalogs == len(catalog_names))]
    return flagged.sort_values("n_members", ascending=False)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--catalogs", nargs="+", required=True,
        help="NAME=path.bed pairs, e.g. GENCODE47=gencode47.exons.bed"
    )
    ap.add_argument("--min-frac", type=float, default=0.3,
                     help="Min fraction of the smaller gene's exonic bp that must be covered (default 0.3)")
    ap.add_argument("--outdir", type=Path, default=Path("consensus_out"))
    ap.add_argument("--dump-components", action="store_true",
                     help="Also write a long-format (locus, catalog, gene_id) table for full manual audit")
    ap.add_argument("--no-dual-bar-plot", action="store_true",
                     help="Skip generating raw_vs_loci_dual_bar.png (requires matplotlib)")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    tmpdir = Path(tempfile.mkdtemp())

    catalogs = {}
    for spec in args.catalogs:
        name, path = spec.split("=", 1)
        catalogs[name] = Path(path)

    print(f"[1/4] Building per-gene exon footprints for {len(catalogs)} catalogs...", file=sys.stderr)
    footprints = {}
    lengths = {}
    footprint_bed_paths = {}
    catalog_gene_ids = {}
    for name, path in catalogs.items():
        fp = build_exon_footprints(path)
        footprints[name] = fp
        lengths[name] = footprint_lengths(fp)
        catalog_gene_ids[name] = set(fp["gene_id"].unique())
        bed_out = tmpdir / f"{name}.footprint.bed"
        write_bed(fp, bed_out)
        footprint_bed_paths[name] = bed_out
        print(f"  {name}: {len(catalog_gene_ids[name])} genes", file=sys.stderr)

    print(f"[2/4] Pairwise bedtools intersect + overlap-fraction filtering (min_frac={args.min_frac})...",
          file=sys.stderr)
    edges_by_pair = {}
    for cat_a, cat_b in itertools.combinations(catalogs.keys(), 2):
        overlap_df = pairwise_overlap_bp(footprint_bed_paths[cat_a], footprint_bed_paths[cat_b])
        if overlap_df.empty:
            edges_by_pair[(cat_a, cat_b)] = overlap_df.assign(frac_smaller=[])
            continue
        edges = build_edges(overlap_df, lengths[cat_a], lengths[cat_b], args.min_frac)
        edges_by_pair[(cat_a, cat_b)] = edges
        print(f"  {cat_a} vs {cat_b}: {len(edges)} gene-gene edges pass threshold", file=sys.stderr)

    print("[3/4] Building overlap graph and taking connected components...", file=sys.stderr)
    matrix_df, components, node_catalog, node_gene = build_consensus_matrix(catalog_gene_ids, edges_by_pair)
    print(f"  {len(matrix_df)} consensus loci", file=sys.stderr)

    matrix_df = add_relationship_type(matrix_df, list(catalogs.keys()))
    matrix_path = args.outdir / "consensus_matrix.tsv"
    matrix_df.to_csv(matrix_path, sep="\t", index=False)
    print(f"  wrote {matrix_path}", file=sys.stderr)
    print("  (columns '<catalog>' = 0/1 presence for UpSet; "
          "'<catalog>_n' = raw gene count absorbed into that locus; "
          "'relationship_type' = 1:1 / 1:many / many:many / singleton per locus)", file=sys.stderr)

    print("[4/4] Flagging likely chaining artifacts (top 1% by component size)...", file=sys.stderr)
    flagged = flag_outlier_components(matrix_df, list(catalogs.keys()))
    flagged_path = args.outdir / "flagged_components_for_review.tsv"
    flagged.to_csv(flagged_path, sep="\t", index=False)
    print(f"  {len(flagged)} components flagged, wrote {flagged_path}", file=sys.stderr)

    split_df = detect_split_events(matrix_df, list(catalogs.keys()))
    split_path = args.outdir / "split_events.tsv"
    split_df.to_csv(split_path, sep="\t", index=False)
    print(f"  {len(split_df)} loci show a >1-gene-from-one-catalog split, wrote {split_path}",
          file=sys.stderr)

    if args.dump_components:
        long_df = dump_components_long(components, node_catalog, node_gene)
        long_path = args.outdir / "components_long_format.tsv"
        long_df.to_csv(long_path, sep="\t", index=False)
        print(f"  wrote full member audit table {long_path}", file=sys.stderr)

    frag_df = fragmentation_index(matrix_df, catalog_gene_ids)
    frag_path = args.outdir / "fragmentation_index.tsv"
    frag_df.to_csv(frag_path, sep="\t", index=False)
    print(f"  wrote {frag_path} "
          "(fragmentation_index > 1 means the UpSet set-size bar UNDERCOUNTS "
          "that catalog's true gene contribution -- see docstring)", file=sys.stderr)

    if not args.no_dual_bar_plot:
        try:
            dual_bar_path = args.outdir / "raw_vs_loci_dual_bar.png"
            plot_dual_bar(frag_df, dual_bar_path)
            print(f"  wrote {dual_bar_path}", file=sys.stderr)
        except ImportError:
            print("  (skipped dual-bar plot: matplotlib not installed -- "
                  "pip install matplotlib --break-system-packages)", file=sys.stderr)

    # sanity check printout -- uses TRUE per-catalog gene counts (sum of *_n
    # columns), not presence sums, so 1-to-many splits no longer look like
    # lost genes here.
    print("\nSanity check -- per-catalog gene counts recovered in the consensus matrix:", file=sys.stderr)
    for name in catalogs:
        n_true = matrix_df[f"{name}_n"].sum()
        n_loci = matrix_df[name].sum()
        raw = len(catalog_gene_ids[name])
        flag = "" if n_true == raw else "  <-- MISMATCH, investigate"
        print(f"  {name}: {n_true} genes recovered across {n_loci} loci "
              f"(raw catalog gene count: {raw}){flag}", file=sys.stderr)


if __name__ == "__main__":
    main()
