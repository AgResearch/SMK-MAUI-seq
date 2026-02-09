#!/usr/bin/env python3
# Copyright: CC BY-SA 4.0 - 2025 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

"""
merge_clr_transform.py

Merge MAUI-seq accepted_by_sample_sequences.tab tables across genes,
convert counts to within-gene proportions (closure), apply CLR
transformation within samples, and write out summary metrics.

Workflow:
    1. Read accepted_by_sample_sequences.tab for each gene
    2. Filter out metadata rows (total, seconds)
    3. Filter out samples with all-zero counts across all genes
    4. Convert counts to within-gene proportions (closure)
    5. Merge gene tables by sample (gene-prefixed column names)
    6. Apply multiplicative zero replacement
    7. CLR transform proportions within each sample (row)
    8. Write merged counts, proportions, CLR-transformed table, and summary
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
        description="Merge MAUI-seq gene count tables and CLR transform.",
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
        help="Directory to write output files.",
    )
    parser.add_argument(
        "--genes",
        nargs="+",
        default=GENES,
        help="Gene amplicons to process.",
    )
    parser.add_argument(
        "--input-table",
        default="accepted_by_sample_sequences.tab",
        help="Name of the per-sample count table inside MAUIcount_output/.",
    )
    parser.add_argument(
        "--pseudocount",
        type=float,
        default=0.5,
        help=(
            "Pseudocount for multiplicative zero replacement prior to CLR. "
            "Applied as a fraction of the minimum non-zero proportion."
        ),
    )
    parser.add_argument(
        "--min-total-counts",
        type=int,
        default=0,
        help="Drop samples with fewer total counts (across all genes) than this.",
    )
    return parser.parse_args(argv)


def read_gene_table(filepath: str) -> pd.DataFrame:
    """Read a MAUIcount accepted_by_sample_sequences.tab file."""
    df = pd.read_csv(filepath, sep="\t", index_col=0)
    # Drop unnamed trailing columns (artefact of trailing tabs)
    df = df.loc[:, ~df.columns.str.contains(r"^Unnamed")]
    # Drop metadata rows
    df = df.drop(index=[r for r in METADATA_ROWS if r in df.index], errors="ignore")
    # Ensure integer counts
    df = df.fillna(0).astype(int)
    return df


def counts_to_proportions(df: pd.DataFrame) -> pd.DataFrame:
    """Convert counts to within-gene proportions (closure) per sample."""
    row_sums = df.sum(axis=1)
    # Avoid division by zero for all-zero samples
    row_sums = row_sums.replace(0, np.nan)
    proportions = df.div(row_sums, axis=0)
    return proportions.fillna(0.0)


def multiplicative_zero_replacement(
    df: pd.DataFrame, pseudocount_fraction: float = 0.5
) -> pd.DataFrame:
    """
    Multiplicative zero replacement (Martín-Fernández et al., 2003).

    For each row, zeros are replaced with delta = pseudocount_fraction *
    (minimum non-zero value in that row). Non-zero values are adjusted
    downward so the row still sums to 1 (closure preserved).

    Delta is capped at 1/D (where D = total number of parts) so that each
    zero receives at most an equal share of the total probability budget,
    ensuring the adjustment factor for non-zero values remains positive.
    """
    result = df.copy()
    D = df.shape[1]  # total number of parts (columns)
    for idx in result.index:
        row = result.loc[idx].values.astype(float)
        nonzero = row[row > 0]
        if len(nonzero) == 0:
            continue  # skip all-zero rows
        n_zeros = (row == 0).sum()
        if n_zeros == 0:
            continue  # no zeros to replace
        delta = pseudocount_fraction * nonzero.min()
        # Cap delta at 1/D: each zero gets at most an equal share of budget
        max_delta = 1.0 / D
        if delta > max_delta:
            delta = max_delta
        # Replace zeros with delta
        zero_mask = row == 0
        nonzero_mask = ~zero_mask
        row[zero_mask] = delta
        # Scale non-zero values so row still sums to 1
        adjustment = 1 - n_zeros * delta
        row[nonzero_mask] = row[nonzero_mask] * adjustment / row[nonzero_mask].sum()
        result.loc[idx] = row
    return result


def clr_transform(df: pd.DataFrame) -> pd.DataFrame:
    """
    Centered log-ratio transformation within each sample (row).

    CLR(x_i) = ln(x_i / g(x)) where g(x) is the geometric mean of the row.
    """
    log_data = np.log(df.values.astype(float))
    geom_mean = log_data.mean(axis=1, keepdims=True)  # per row
    clr_data = log_data - geom_mean
    return pd.DataFrame(clr_data, index=df.index, columns=df.columns)


def compute_summary(
    merged_counts: pd.DataFrame,
    merged_proportions: pd.DataFrame,
    clr_df: pd.DataFrame,
    gene_tables: dict,
) -> str:
    """Generate a human-readable summary of the merge and transformation."""
    lines = [
        "MAUI-seq CLR Transformation Summary",
        "=" * 50,
        "",
        f"Genes processed: {', '.join(gene_tables.keys())}",
        f"Samples: {len(merged_counts)}",
        f"Total seq types (parts): {len(merged_counts.columns)}",
        "",
    ]

    for gene, df in gene_tables.items():
        gene_cols = [c for c in merged_counts.columns if c.startswith(f"{gene}_")]
        total = merged_counts[gene_cols].sum().sum()
        lines.append(f"  {gene}: {len(gene_cols)} seq types, {total:,} total counts")

    lines.append("")
    lines.append("Per-sample total counts (across all genes):")
    row_totals = merged_counts.sum(axis=1).sort_values(ascending=False)
    for sample, total in row_totals.items():
        lines.append(f"  {sample}: {total:,}")

    lines.append("")
    lines.append("Zero counts in merged proportions table:")
    total_cells = merged_proportions.size
    zero_cells = (merged_proportions == 0).sum().sum()
    lines.append(f"  {zero_cells} / {total_cells} ({100 * zero_cells / total_cells:.1f}%)")

    lines.append("")
    lines.append("CLR value range:")
    lines.append(f"  Min: {clr_df.min().min():.4f}")
    lines.append(f"  Max: {clr_df.max().max():.4f}")
    lines.append(f"  Mean: {clr_df.values.mean():.6f} (should be ~0)")

    # Samples that were all-zero in any gene
    lines.append("")
    lines.append("Samples with zero counts per gene:")
    for gene, df in gene_tables.items():
        zero_samples = df.index[df.sum(axis=1) == 0].tolist()
        if zero_samples:
            lines.append(f"  {gene}: {', '.join(zero_samples)}")

    return "\n".join(lines)


def main(argv=None):
    args = parse_args(argv)

    os.makedirs(args.output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Read gene tables
    # ------------------------------------------------------------------
    gene_tables = {}
    for gene in args.genes:
        filepath = os.path.join(
            args.input_dir, gene, "MAUIcount_output", args.input_table
        )
        if not os.path.isfile(filepath):
            # Fall back to accepted_sequences.tab
            alt_filepath = os.path.join(
                args.input_dir, gene, "MAUIcount_output", "accepted_sequences.tab"
            )
            if os.path.isfile(alt_filepath):
                logger.warning(
                    "Could not find %s for %s, using accepted_sequences.tab",
                    args.input_table,
                    gene,
                )
                filepath = alt_filepath
            else:
                logger.error(
                    "No count table found for gene %s in %s", gene, args.input_dir
                )
                sys.exit(1)

        logger.info("Reading %s", filepath)
        gene_tables[gene] = read_gene_table(filepath)

    # ------------------------------------------------------------------
    # 2. Prefix columns and merge by sample
    # ------------------------------------------------------------------
    frames_counts = []
    frames_proportions = []

    for gene, df in gene_tables.items():
        prefixed = df.copy()
        prefixed.columns = [f"{gene}_{col}" for col in prefixed.columns]
        frames_counts.append(prefixed)

        props = counts_to_proportions(df)
        props.columns = [f"{gene}_{col}" for col in props.columns]
        frames_proportions.append(props)

    # Outer join so samples missing from one gene get zeros
    merged_counts = pd.concat(frames_counts, axis=1).fillna(0).astype(int)
    merged_proportions = pd.concat(frames_proportions, axis=1).fillna(0.0)

    # Sort samples alphabetically
    merged_counts = merged_counts.sort_index()
    merged_proportions = merged_proportions.sort_index()

    # ------------------------------------------------------------------
    # 3. Filter samples
    # ------------------------------------------------------------------
    # Drop samples where ALL counts are zero
    all_zero_mask = merged_counts.sum(axis=1) == 0
    if all_zero_mask.any():
        dropped = merged_counts.index[all_zero_mask].tolist()
        logger.warning("Dropping all-zero samples: %s", dropped)
        merged_counts = merged_counts[~all_zero_mask]
        merged_proportions = merged_proportions[~all_zero_mask]

    # Drop samples below minimum total count threshold
    if args.min_total_counts > 0:
        low_mask = merged_counts.sum(axis=1) < args.min_total_counts
        if low_mask.any():
            dropped = merged_counts.index[low_mask].tolist()
            logger.warning(
                "Dropping samples below %d total counts: %s",
                args.min_total_counts,
                dropped,
            )
            merged_counts = merged_counts[~low_mask]
            merged_proportions = merged_proportions[~low_mask]

    logger.info(
        "Merged table: %d samples x %d seq types",
        merged_counts.shape[0],
        merged_counts.shape[1],
    )

    # ------------------------------------------------------------------
    # 4. Zero replacement and CLR transformation
    # ------------------------------------------------------------------
    logger.info("Applying multiplicative zero replacement (fraction=%.2f)", args.pseudocount)
    replaced = multiplicative_zero_replacement(merged_proportions, args.pseudocount)

    logger.info("Applying CLR transformation within samples (across columns)")
    clr_df = clr_transform(replaced)

    # ------------------------------------------------------------------
    # 5. Write outputs
    # ------------------------------------------------------------------
    counts_path = os.path.join(args.output_dir, "merged_counts.tsv")
    proportions_path = os.path.join(args.output_dir, "merged_proportions.tsv")
    replaced_path = os.path.join(args.output_dir, "merged_proportions_zero_replaced.tsv")
    clr_path = os.path.join(args.output_dir, "merged_clr.tsv")
    summary_path = os.path.join(args.output_dir, "clr_summary.txt")

    merged_counts.to_csv(counts_path, sep="\t")
    logger.info("Wrote %s", counts_path)

    merged_proportions.to_csv(proportions_path, sep="\t", float_format="%.6f")
    logger.info("Wrote %s", proportions_path)

    replaced.to_csv(replaced_path, sep="\t", float_format="%.6f")
    logger.info("Wrote %s", replaced_path)

    clr_df.to_csv(clr_path, sep="\t", float_format="%.6f")
    logger.info("Wrote %s", clr_path)

    summary = compute_summary(merged_counts, merged_proportions, clr_df, gene_tables)
    with open(summary_path, "w") as f:
        f.write(summary + "\n")
    logger.info("Wrote %s", summary_path)

    print(summary)

    return 0


if __name__ == "__main__":
    sys.exit(main())
