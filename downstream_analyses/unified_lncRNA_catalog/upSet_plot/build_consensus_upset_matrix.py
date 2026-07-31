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
def build_consensus_matrix(catalog_gene_ids: dict, edges_by_pair: dict) -> pd.DataFrame:
    """
    catalog_gene_ids: {catalog_name: set(gene_ids)}
    edges_by_pair:    {(cat_a, cat_b): edge_df with columns gene_a, gene_b}
    """
    G = nx.Graph()

    # every gene is a node, tagged with its catalog of origin
    node_catalog = {}
    for cat, genes in catalog_gene_ids.items():
        for g in genes:
            node_id = f"{cat}::{g}"  # namespace gene ids by catalog to avoid collisions
            G.add_node(node_id)
            node_catalog[node_id] = cat

    for (cat_a, cat_b), edge_df in edges_by_pair.items():
        for _, row in edge_df.iterrows():
            G.add_edge(f"{cat_a}::{row['gene_a']}", f"{cat_b}::{row['gene_b']}")

    components = list(nx.connected_components(G))

    matrix_rows = []
    for i, comp in enumerate(components):
        locus_id = f"locus_{i:06d}"
        present = {cat: 0 for cat in catalog_gene_ids}
        for node in comp:
            present[node_catalog[node]] = 1
        row = {"consensus_locus": locus_id, "n_members": len(comp)}
        row.update(present)
        matrix_rows.append(row)

    return pd.DataFrame(matrix_rows)


# ---------------------------------------------------------------------------
# Diagnostics: flag likely chaining artifacts
# ---------------------------------------------------------------------------
def flag_outlier_components(matrix_df: pd.DataFrame, catalog_names: list, pct=99):
    """
    Components with unusually many member genes or unusually many catalogs
    are worth manual review -- likely bridge/chaining artifacts even after
    switching to exon-footprint overlap.
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
    matrix_df = build_consensus_matrix(catalog_gene_ids, edges_by_pair)
    print(f"  {len(matrix_df)} consensus loci", file=sys.stderr)

    matrix_path = args.outdir / "consensus_matrix.tsv"
    matrix_df.to_csv(matrix_path, sep="\t", index=False)
    print(f"  wrote {matrix_path}", file=sys.stderr)

    print("[4/4] Flagging likely chaining artifacts (top 1% by component size)...", file=sys.stderr)
    flagged = flag_outlier_components(matrix_df, list(catalogs.keys()))
    flagged_path = args.outdir / "flagged_components_for_review.tsv"
    flagged.to_csv(flagged_path, sep="\t", index=False)
    print(f"  {len(flagged)} components flagged, wrote {flagged_path}", file=sys.stderr)

    # sanity check printout
    print("\nSanity check -- per-catalog set sizes in the consensus matrix:", file=sys.stderr)
    for name in catalogs:
        n = matrix_df[name].sum()
        print(f"  {name}: {n} (raw catalog gene count: {len(catalog_gene_ids[name])})", file=sys.stderr)


if __name__ == "__main__":
    main()
