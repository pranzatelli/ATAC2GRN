# HMMtool

HMMtool segments Tn5-insertion signal from ATAC2GRN `Pipes` output into
discrete chromatin states using a Poisson-emission Hidden Markov Model, and
optionally correlates each state's mean signal against a gene's RNA-seq
expression across samples.

This is a standalone downstream analysis script, structured like `ROCtool/`
(a documented Python script, run directly — see "Why no wrapper script?"
below for why it has no `.sh` companion).

## What it does

1. For a chosen genomic region (`--chrom`/`--start`/`--end`), builds a
   per-base Tn5-insertion count track for every sample/project under an
   ATAC2GRN `Pipes` run's `Output/` directory (using each project's
   `{project}-picard.bam`, with the same Tn5 +4/-5 offset correction used
   elsewhere in ATAC2GRN).
2. Optionally (see below) correlates a rolling-window-smoothed version of
   that signal against a gene/transcript's RNA-seq expression across the
   samples that have both ATAC and RNA-seq data.
3. Fits a `PoissonHMM` (an `hmmlearn`-based Poisson-emission HMM) to the
   smoothed multi-sample signal and predicts a discrete state per bin.
4. Writes the predicted states as a bedgraph, plus `.npy` dumps of the raw
   per-bin state calls, fitted Poisson rate parameters, and transition
   matrix.

## Required inputs

- **Must be run after an ATAC2GRN `Pipes` run** (see `../Pipes/README.md`).
  `--atac-output-dir` should point at that run's `Output/` directory, which
  is expected to contain one subfolder per sample/project, each holding a
  `{project}-picard.bam` (and its index).
- A genomic region: `--chrom`, `--start`, `--end` (0-based, half-open, like
  a BED interval).
- A `--label` for this run; outputs are written to `{output-dir}/{label}/`.

## Optional: RNA-seq expression correlation

Expression correlation is entirely opt-in. Provide **all three** of the
following to enable it, or omit all three to skip it and just do ATAC-only
segmentation:

- `--genes` — comma-separated gene/transcript ID(s) (matching the IDs used
  in the RNA-seq quantification files, e.g. `quant.sf` transcript IDs;
  a trailing version suffix like `.5` is stripped automatically). Their
  TPM values are summed per sample before correlation.
- `--expression-quant-dir` — a directory with one subfolder per RNA-seq
  sample, each containing a [salmon](https://combine-lab.github.io/salmon/)
  `quant.sf` file.
- `--sample-metadata` — a TSV/CSV file mapping RNA-seq sample IDs to
  ATAC-seq sample/project IDs, since the two are rarely named identically.
  By default it must have two columns named `rna_sample` and `atac_sample`;
  override the column names with `--metadata-rna-col`/`--metadata-atac-col`
  if your file uses different headers. Example:

  ```csv
  rna_sample,atac_sample
  RNA_Sample_01,ATAC_Sample_01
  RNA_Sample_02,ATAC_Sample_02
  ```

  Only samples present in this mapping (and with both an ATAC BAM and an
  RNA-seq `quant.sf`) are used for correlation.

This generalizes the original script's ENCODE-specific metadata lookup
(which joined RNA-seq and ATAC-seq experiment accessions indirectly through
ENCODE's `Biosample term name`/`Age`/`Assay` columns) into a direct,
institution-agnostic sample-to-sample mapping. If you need to derive such a
mapping from ENCODE metadata yourself, build the two-column CSV above from
it first.

## CLI usage

```
python HMMtool.py \
  --chrom chr11 --start 47168187 --end 51230222 \
  --label FOLH1_region \
  --atac-output-dir /path/to/Pipes/Output \
  --n-states 25 --n-iter 50 --rolling-window 101
```

With optional RNA-seq correlation:

```
python HMMtool.py \
  --chrom chr11 --start 47168187 --end 51230222 \
  --label FOLH1_region \
  --atac-output-dir /path/to/Pipes/Output \
  --genes NM_001014986,NM_004476 \
  --expression-quant-dir /path/to/salmon_quant \
  --sample-metadata sample_metadata.csv \
  --output-dir Output --n-states 25 --n-iter 50 --rolling-window 101
```

Run `python HMMtool.py --help` for the full flag list, including
`--exclude-samples` (comma-separated sample/project names to drop before
fitting — no names are hardcoded; supply your own if you need to exclude
known-bad samples/runs) and `--random-seed`.

## Output files

Written to `{output-dir}/{label}/` (default `output-dir` is `Output`):

| File | Always written? | Contents |
|---|---|---|
| `{label}.bedgraph` | yes | Predicted chromatin state per genomic interval |
| `scores.npy` | yes | Per-bin predicted state (same order as the bedgraph) |
| `lambdas.npy` | yes | Fitted per-state, per-sample Poisson rate parameters |
| `transmat.npy` | yes | Fitted HMM state-transition matrix |
| `corr.npy` | only with expression correlation | Per-bin correlation of smoothed ATAC signal with gene expression |
| `{label}.csv` | only with expression correlation | Summed per-sample expression for the requested gene(s) |

## Why no wrapper script?

`ROCtool.sh` exists because ROCtool needs a `bedtools`-based preprocessing
step before invoking Python. HMMtool doesn't need any shell preprocessing —
all inputs are directories/files taken directly as CLI arguments — so this
folder intentionally has no `HMMtool.sh`. Flagging this as a deliberate
deviation from the `ROCtool/` convention rather than an oversight.

## Notes on the refactor

This tool is a refactor of a user-supplied research script. Besides adding
the CLI and generalizing the RNA-seq metadata handling described above, two
bugs were fixed along the way:

- `_generate_sample_from_state` referenced an undefined variable; it now
  correctly samples from that state's fitted Poisson rates.
- The bedgraph was previously written using bin indices from *before*
  the rolling-window smoothing step, while the HMM states are predicted
  *after* smoothing (which can drop a few edge bins) — a latent
  off-by-a-few-bins mismatch. The bedgraph now uses the indices the states
  are actually aligned to.

The core `PoissonHMM` model and the Tn5-offset BAM-to-array extraction
(`construct_atac_array`) are unchanged from the original, validated
implementation.

`hmmlearn` is pinned to `0.3.3` in `../atac2grn.def` because `PoissonHMM`
subclasses `hmmlearn`'s private `_BaseHMM` API, whose method signatures are
not guaranteed stable across releases; this is the version HMMtool's
overrides were verified against.
