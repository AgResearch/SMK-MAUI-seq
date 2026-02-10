#!/usr/bin/env python3
# Copyright: CC BY-SA 4.0 - 2025 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

"""
compose_sample_fasta.py

For each gene, compose a FASTA file containing one entry per sample per
accepted sequence type that is present (count > 0) in that sample.

Each FASTA header includes the sample name, gene, seq type identifier,
raw count, and within-sample proportion.

Input:
    - accepted_sequences.fas      (seq_id -> nucleotide sequence)
    - accepted_by_sample_sequences.tab  (sample x seq_id count matrix)

Output:
    - One FASTA file per gene: <gene>_by_sample.fasta
"""

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

METADATA_ROWS = {"total", "seconds"}
GENES = ["recA", "rpoB", "nodA", "nodD"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Compose per-sample FASTA files from MAUIcount accepted "
            "sequences and per-sample count tables."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help=(
            "Path to the MAUIcount output directory containing gene "
            "subdirectories (e.g., results/<RUN>/04_MAUIcount)."
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to write output FASTA files.",
    )
    parser.add_argument(
        "--genes",
        nargs="+",
        default=GENES,
        help="Gene amplicons to process.",
    )
    parser.add_argument(
        "--count-table",
        default="accepted_by_sample_sequences.tab",
        help="Name of the per-sample count table inside MAUIcount_output/.",
    )
    parser.add_argument(
        "--fasta-file",
        default="accepted_sequences.fas",
        help="Name of the accepted sequences FASTA inside MAUIcount_output/.",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="Minimum count to include a sequence for a sample.",
    )
    return parser.parse_args(argv)


def read_fasta(filepath: str) -> dict:
    """
    Read a FASTA file and return a dict of {seq_id: sequence}.

    Headers are expected in the format: >seq_N_totalcount_seconds
    The seq_id extracted is 'seq_N' (first two underscore-delimited fields).
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Write previous entry
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                # Parse header: >seq_1_22240_577 -> seq_1
                header = line[1:]  # strip '>'
                parts = header.split("_")
                current_id = f"{parts[0]}_{parts[1]}"
                current_seq = []
            else:
                current_seq.append(line)

    # Last entry
    if current_id is not None:
        sequences[current_id] = "".join(current_seq)

    return sequences


def read_count_table(filepath: str) -> pd.DataFrame:
    """Read a MAUIcount accepted_by_sample_sequences.tab file."""
    df = pd.read_csv(filepath, sep="\t", index_col=0)
    # Drop unnamed trailing columns (artefact of trailing tabs)
    df = df.loc[:, ~df.columns.str.contains(r"^Unnamed")]
    # Drop metadata rows
    df = df.drop(index=[r for r in METADATA_ROWS if r in df.index], errors="ignore")
    # Ensure integer counts
    df = df.fillna(0).astype(int)
    return df


def compose_gene_fasta(
    gene: str,
    sequences: dict,
    count_table: pd.DataFrame,
    min_count: int = 1,
) -> list:
    """
    Compose FASTA entries for a single gene.

    Returns a list of (header, sequence) tuples.
    """
    entries = []

    for sample in sorted(count_table.index):
        row = count_table.loc[sample]
        row_total = row.sum()
        if row_total == 0:
            continue  # skip all-zero samples

        for seq_id in count_table.columns:
            count = int(row[seq_id])
            if count < min_count:
                continue

            proportion = count / row_total
            seq = sequences.get(seq_id)
            if seq is None:
                logger.warning(
                    "Sequence %s not found in FASTA for gene %s", seq_id, gene
                )
                continue

            header = (
                f">gene={gene};seq={seq_id};sample={sample};"
                f"count={count};proportion={proportion:.6f}"
            )
            entries.append((header, seq))

    return entries


def main(argv=None):
    args = parse_args(argv)

    os.makedirs(args.output_dir, exist_ok=True)

    total_entries = 0

    for gene in args.genes:
        fasta_path = os.path.join(
            args.input_dir, gene, "MAUIcount_output", args.fasta_file
        )
        table_path = os.path.join(
            args.input_dir, gene, "MAUIcount_output", args.count_table
        )

        if not os.path.isfile(fasta_path):
            logger.error("FASTA not found: %s", fasta_path)
            sys.exit(1)
        if not os.path.isfile(table_path):
            logger.error("Count table not found: %s", table_path)
            sys.exit(1)

        logger.info("Reading %s", fasta_path)
        sequences = read_fasta(fasta_path)
        logger.info("  %d sequences loaded", len(sequences))

        logger.info("Reading %s", table_path)
        count_table = read_count_table(table_path)
        logger.info(
            "  %d samples x %d seq types",
            count_table.shape[0],
            count_table.shape[1],
        )

        entries = compose_gene_fasta(gene, sequences, count_table, args.min_count)
        logger.info("  %d FASTA entries composed", len(entries))

        output_path = os.path.join(args.output_dir, f"{gene}_by_sample.fasta")
        with open(output_path, "w") as f:
            for header, seq in entries:
                f.write(f"{header}\n{seq}\n")

        logger.info("Wrote %s", output_path)
        total_entries += len(entries)

    logger.info(
        "Done. %d total entries across %d genes.", total_entries, len(args.genes)
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
