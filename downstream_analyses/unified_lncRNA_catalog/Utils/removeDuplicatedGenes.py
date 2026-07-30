#!/usr/bin/env python3

"""
Written to get deduplicated (not intersecting with each other)
list of geneIDs from raw bedtools intersect output:
bedtools intersect -s -wa -wb -f 1.00 -F 1.00 -e -a file1.bed -b file1.bed
e.g.: chr2	220791077	220797893	bigTranscriptome|LOC_000000005443	0	+	chr2	220791077	220797893	gen27|ENSG00000224819.1	0	+
We aim to select 1 geneID from group of mutually intersecting genes.
Obviously, output can contain entries where geneID1 is intersecting the same geneID1.
Input: bedtools output as file or STDIN
Output: list of genes in STDOUT
"""

import sys
import argparse
from contextlib import nullcontext


def main() -> None:
    """
    Handles input from STDIN.
    Returns geneIDs to STDOUT.
    """
    # Handle comandline arguments
    parser = argparse.ArgumentParser(
        description="Extracts deduplicated geneIDs from bedtools intersect output with mutually intersecting genes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage="\n%(prog)s input.bed\ncat input.tsv | %(prog)s\n%(prog)s < input.tsv",
    )
    parser.add_argument(
        "input", nargs="?", default="-", help="Input file ('-' or omitted means stdin)"
    )

    args = parser.parse_args()

    # Check if input is a file or STDIN
    if args.input == "-" and sys.stdin.isatty():
        parser.error("No input provided. Specify an input file or pipe data to stdin.")
    elif args.input == "-":
        # Trick to not close STDIN on completion
        cm = nullcontext(sys.stdin)
    else:
        # If it's file - it will be closed on completion
        cm = open(args.input, "r")

    # Create sets
    kept_genes = set()
    duplicates = set()
    possibly_specific = set()

    # Process intput lines
    with cm as infile:
        for i, line in enumerate(infile, start=1):
            # Process line
            features_list = line.rstrip("\n").split(sep="\t")
            # Verify line length
            if len(features_list) < 10:
                parser.error(f"Input line is too short!\nLine: {features_list}")
            # If seems correct, extract genes
            gene1 = features_list[3]
            gene2 = features_list[9]

            # Keep IDs of genes intersecting themselves for future
            if gene1 == gene2:
                possibly_specific.add(gene1)
                continue

            # If geneID is not the duplicate; add to kept_genes
            if gene1 not in duplicates:
                kept_genes.add(gene1)
                duplicates.add(gene2)

            # Not all duplicates can align with each other
            # Add genes2 to duplicates if gene1 is there
            elif (gene1 in duplicates) and (gene2 not in kept_genes):
                duplicates.add(gene2)

            # Do not add kept genes to duplicates
            else:
                pass

            # Report progress
            if i % 1000 == 0:
                print(f"\rProcessed {i:,} lines.", end="", file=sys.stderr, flush=True)
    # Update progress, so ppl won't think last lines were skipped
    print(f"\rProcessed {i:,} lines.", end="", file=sys.stderr, flush=True)

    # Remove genes known as duplicates from possibly catalog specific geneIDs
    possibly_specific.difference_update(duplicates)

    # Add catalog sepcific genes with deduplicated geneIDs
    kept_genes.update(possibly_specific)

    # Write deduplicated genes to STDOUT
    sys.stdout.write("\n".join(kept_genes))
    sys.stdout.write("\n")

    # Report completeness
    print("\nDone!", file=sys.stderr)


if __name__ == "__main__":
    main()
