
# **Fork of MAUI-seq for use at AgResearch**

### MAUIsortgenes.py - deprecated when run as a Snakemake workflow. Cutadapt is used in place to demultiplex the gene specific amplicon sequencing reads after paired-read assembly using PEAR.
### MC_parameters.py - deprecated. Amplicon specific parameters are now passed via an amplicon_config.yaml file.
### MAUIcount.py - refactored as a command-line utility for integration with Snakemake. The main analysis code has been left untouched and encapsulated in main().

```
usage: python MAUIcount.py [-h] [--working-folder WORKING_FOLDER] [--amplicon-config AMPLICON_CONFIG] [--fastq-file-list FASTQ_FILE_LIST] [--read-diff READ_DIFF]
                    [--reject-threshold REJECT_THRESHOLD] [--add-limit ADD_LIMIT] [--output-types {3,4}] [--output-read-counts] --gene
                    {recA,rpoB,nodA,nodD,viciae_nodD}

MAUIcount: Analyze amplicon sequences with UMI-based error correction

This script reads amplicon sequences from a set of fastq files.
The first UMI_len bases of each read are a random tag (UMI or Unique Molecular Identifier).
The script keeps track of how many times each UMI is used with each unique sequence.
For the set of samples, outputs are files with a list of fasta sequences in descending rank
order of abundance, and corresponding tables with the counts of each sequence in each sample.

Three output sets are produced by default:

"accepted_sequences": the main MAUI-seq output, using secondary sequences to filter out errors

"all_primary_sequences": the same UMI counts before error filtering

"read_sequences": 'conventional' analysis of the sequences, ignoring UMIs

An additional output set can be produced by setting output_types = 4. This has the same 
filtering as accepted_sequences, but applied on a per-sample basis, rather than on totals 
across all samples. If allele frequencies vary greatly across samples, this would be 
preferable in principle, but it is not reliable unless read counts are very high. Otherwise,
it can lead to sequences being stochastically deleted from some samples but not others.

In all cases, the outputs are truncated to discard rare sequences that would have
frequencies below add_limit in the overall set of samples. 

Additional output files are summary.txt (various data about the run) and UMI_stats.tab 
(information on the distribution of read numbers per UMI that may be useful for optimising
the protocol).

A file amplicon_config.yaml specifies the files to be analysed and gene-specific parameters.

If there is a file fastq_file_list.txt in the same folder as the data, only the files listed 
in this file will be included in the analysis. If this file is not present, it will be  
created with a list of all files that have the extension .fastq.

options:
  -h, --help            show this help message and exit
  --working-folder WORKING_FOLDER
                        Path to folder containing FASTQ files.
  --amplicon-config AMPLICON_CONFIG
                        Path to amplicon configuration file.
  --fastq-file-list FASTQ_FILE_LIST
                        Path to a file listing FASTQ files to process. If not provided, all .fastq files in the working folder will be used.
  --read-diff READ_DIFF
                        Count UMI only if most abundant sequence has at least read_diff more reads than the next.
  --reject-threshold REJECT_THRESHOLD
                        Reject sequences that occur as second sequences with UMIs at least reject_threshold times as often as they occur as the primary
                        sequence.
  --add-limit ADD_LIMIT
                        Sequences are included, in rank order, until the next would add a fraction less than add_limit (i.e. this discards sequences with
                        an overall relative abundance less than add_limit).
  --output-types {3,4}  If set to 4, additional output is produced based on filtering separately for each sample.
  --output-read-counts  Output the total number of reads contributing to the UMI primary and secondary counts.
  --gene {recA,rpoB,nodA,nodD,viciae_nodD}
                        Gene amplicon to analyze.

Written by Peter Young. Version 1.01 on 20 January 2020.
Refactored by Ben Perry. Version 2.0 on 10 June 2025.
```

A conda environment config for MAUIcount.py can be found in `workflow/envs/mauiseq.yaml`. Amplicon configuration parameters must be defined in a yaml file, and passed via the CLI. The production config used by the Snakemake workflow is `resources/amplicon_config_production.yaml`; a test config is at `resources/amplicon_config_test.yaml`. Configuration must appear:

```
genes:
  recA:
    UMI_len: 12
    f_primer_len: 26
    r_primer_len: 23
    total_len: 311
  
[...]
```

**Note:** The Snakemake workflow runs MAUIcount with `--output-types 4` so that per-sample count tables (`accepted_by_sample_sequences.tab`) are produced for downstream CLR transformation.

Test data and results have been moved to `resources/Test` and are still valid to check installation and usage. Test data can be used with the following command after activating the mauiseq.yaml conda environment:

```
python workflow/scripts/MAUIcount.py --working-folder resources/Test --amplicon-config resources/amplicon_config_test.yaml --gene recA
```

### Snakemake Workflow

The workflow is defined in `workflow/Snakefile` and requires Snakemake ≥ 9. It is configured via `config/pipeline_config.yaml` (set `RUN_PATH`) or by passing `--config RUN_PATH=<path>` on the command line. Example:

```
snakemake --profile config/slurm_profiles/eri \
  --config RUN_PATH=/path/to/sequencing/run \
  --snakefile workflow/Snakefile
```

The workflow produces the following output structure under `results/<RUN>/`:

| Directory | Rule(s) | Contents |
|---|---|---|
| `00_QC/` | `report_seqkit_*` | SeqKit statistics for raw, assembled, discarded, and unassembled reads |
| `01_assembled/` | `pear_read_assemble` | PEAR paired-end assembled FASTQ files |
| `02_cutadapt/` | `cutadapt_amplicon_demux` | Cutadapt demultiplexed reads per sample |
| `03_amplicons/` | `gather_gene_amplicons` | Per-gene FASTQ files (recA, rpoB, nodA, nodD) |
| `04_MAUIcount/` | `run_MAUIcount` | MAUIcount UMI-deduplicated count tables per gene |
| `05_CLR/` | `merge_clr_transform` | Merged compositions and CLR-transformed output |

Update B.J. Perry February 2026

### merge_clr_transform.py — Compositional Merge and CLR Transformation

This script merges the per-sample UMI-deduplicated count tables from all four gene amplicons (recA, rpoB, nodA, nodD) into a single compositional dataset and applies a centered log-ratio (CLR) transformation for downstream statistical analysis.

#### Processing Steps

1. **Read** — Loads `accepted_by_sample_sequences.tab` from each gene's `MAUIcount_output/` directory. Each table has samples as rows and seq types as columns, with raw UMI-corrected counts.

2. **Filter metadata rows** — Drops the `total` and `seconds` summary rows appended by MAUIcount.

3. **Closure (within-gene)** — Converts each gene's count table to proportions independently per sample (each row sums to 1 within that gene). This removes between-gene sequencing depth bias so that a gene with higher total reads does not dominate the composition.

4. **Merge by sample** — Concatenates the four proportion tables column-wise with gene-prefixed names (e.g. `recA_seq_1`, `rpoB_seq_1`, …), producing a single samples × seq-types table. Samples missing from a gene receive zeros.

5. **Filter samples** — Drops samples that are all-zero across every gene, and optionally those below a `--min-total-counts` threshold.

6. **Multiplicative zero replacement** ([Martín-Fernández et al., 2003](https://doi.org/10.1023/A:1023866030544)) — For each sample (row), zeros are replaced with δ = 0.5 × min(non-zero values in that row), capped at 1/D where D = total parts. Non-zero values are scaled down proportionally so the row still sums to 1. This is necessary because CLR requires strictly positive values.

7. **CLR transformation** — Within each sample (row), computes CLR(x_i) = ln(x_i) − mean(ln(x)) where the mean is the log of the geometric mean across all parts. Each transformed row sums to exactly 0.

#### Outputs (5 files written to `--output-dir`)

| File | Contents |
|---|---|
| `merged_counts.tsv` | Raw integer counts, merged across genes |
| `merged_proportions.tsv` | Within-gene proportions before zero replacement |
| `merged_proportions_zero_replaced.tsv` | Proportions after multiplicative zero replacement (CLR input) |
| `merged_clr.tsv` | CLR-transformed values (the primary analysis-ready output) |
| `clr_summary.txt` | Sample counts, gene breakdown, sparsity, CLR value range, per-gene zero-sample list |

#### CLI Usage

```
usage: merge_clr_transform.py [-h] --input-dir INPUT_DIR --output-dir OUTPUT_DIR
                              [--genes GENES [GENES ...]]
                              [--input-table INPUT_TABLE]
                              [--pseudocount PSEUDOCOUNT]
                              [--min-total-counts MIN_TOTAL_COUNTS]

  --input-dir       Path to 04_MAUIcount directory (contains gene subdirs)
  --output-dir      Where to write the 5 output files
  --genes           Gene amplicons to process (default: recA rpoB nodA nodD)
  --input-table     Table filename inside MAUIcount_output/ (default: accepted_by_sample_sequences.tab)
  --pseudocount     Zero replacement fraction ε (default: 0.5)
  --min-total-counts  Drop samples below this total count (default: 0)
```

This step is integrated into the Snakemake workflow as rule `merge_clr_transform` and produces outputs in `results/<RUN>/05_CLR/`. It requires `--output-types 4` to be set in the `run_MAUIcount` rule so that `accepted_by_sample_sequences.tab` is generated.

Update B.J. Perry February 2026

---
---
---

MAUI-seq: Metabarcoding using amplicons with unique molecular identifiers to improve error correction
===

This is a method for using PCR amplicons to describe microbial diversity, described in Fields, B, Moeskjær, S, Friman, V‐P, Andersen, SU, Young, JPW. MAUI‐seq: Metabarcoding using amplicons with unique molecular identifiers to improve error correction. Mol Ecol Resour. 2020; 00: 1– 18. https://doi.org/10.1111/1755-0998.13294 

This repository includes scripts used in that paper for the analysis of amplicon sequences that incorporate a random sequence tag (UMI or Unique Molecular Identifier).



**MAUIsortgenes** is a simple script that sorts multiplexed fastq sequences into gene-specific files based on a short motif in the forward primer. It was written for our own data but is easily modified.

**MAUIcount** analyses amplicon data for one gene across multiple samples, producing a table of allele frequencies and a fasta file of allele sequences.

**MC_parameters** is a user-edited file that specifies the data source and parameters for MAUIcount.

Sample data files for testing are provided in the **'Test'** folder.


MAUIsortgenes.py
--- 
This script uses BioPython, but is otherwise written in standard Python 3. The input files contain amplicon sequences in fastq format, for example, those that have been assembled from paired-end reads using PEAR [Zhang, J, Kobert, K, Flouri, T, and Stamatakis, A. (2014). PEAR: a fast and accurate Illumina paired-end read merger. Bioinformatics. 30(5); 614.]. 

It sorts sequences into separate files for each gene, based on a short unambiguous motif (tag) in the forward primer. Sequences that do not match any tag are filed separately as "short" or "unknown". It is very simple and does not allow mismatches. It produces files that are suitable as input for MAUIcount, but users may prefer to use other methods to create these single-gene fastq files.


MAUIcount.py
---
This is the core script for the MAUI method. It is written in standard Python 3 and uses only modules in the Python standard library. It reads amplicon sequences from a set of fastq files.

The first UMI_len bases of each read are a random tag, called a Unique Molecular Identifier (UMI). The script keeps track of how many times each UMI is used with each unique sequence. For the set of samples, outputs are files with a list of fasta sequences in descending rank order of abundance, and corresponding tables with the counts of each sequence in each sample.

Three output sets are produced by default:

"accepted_sequences": the main MAUI-seq output, using secondary sequences to filter out errors

"all_primary_sequences": the same UMI counts before error filtering

"read_sequences": 'conventional' analysis of the sequences, ignoring UMIs

An additional output set can be produced by setting output_types = 4. This has the same filtering as accepted_sequences, but applied on a per-sample basis, rather than on totals across all samples. If allele frequencies vary greatly across samples, this would be preferable in principle, but it is not reliable unless read counts are very high. Otherwise, it can lead to sequences being stochastically deleted from some samples but not others.

In all cases, the outputs are truncated to discard rare sequences that would have frequencies below add_limit in the overall set of samples. 

Additional output files are summary.txt (various data about the run) and UMI_stats.tab (information on the distribution of read numbers per UMI that may be useful for optimising the protocol).

If there is a file fastq_file_list.txt in the same folder as the data, only the files listed in this file will be included in the analysis. If this file is not present, it will be created with a list of all files that have the extension .fastq.


MC_parameters.py
---
This file is imported by MAUIcount.py and specifies the files to be analysed and gene-specific parameters. It should be in the same folder as the MAUIcount.py script (or another path that will be found). It needs to be edited appropriately for your project, and avoids the need for a long list of command line arguments.

Sets of parameters for different genes and data sets can be 'stored' as commented-out lines in MC_parameters.


Test
---
The folder Test provides some sample input data for MAUIcount.py (a subset of the data in Figure 3 of the publication). The supplied version of MC_parameters is set up to run an analysis of these data, provided that MAUIcount is run from the folder immediately above the Test folder. The command "python3 MAUIcount.py" should be sufficient, since all parameters are already set up in MC_parameters.py. If it is working correctly, MAUIcount should generate a new folder called MAUIcount_output within Test that has exactly the same files as the existing folder MAUIcount_expected_output. N.B This software requires Python 3; it has not been tested with Python 2 and may not run as expected.

---
Created by Peter Young peter.young@york.ac.uk


