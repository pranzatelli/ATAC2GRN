# Pipes

This folder contains the pipeline code used in the paper: raw bash scripts
(`bash/`) used during pipeline exploration, and a Snakemake pipeline
(`Snakemake/`) for running a finalized configuration end-to-end. This README
covers the Snakemake pipeline, which is the actively maintained, runnable
entry point for new users.

## Quickstart

### 1. Get the container

The pipeline expects `bowtie2`, `samtools`, `picard.jar`, HOMER
(`makeTagDirectory`/`findPeaks`/`pos2bed.pl`), `bedtools`, `bedops`
(`bedmap`), and `rgt-hint` all on `PATH` (or at documented fixed locations).
The easiest way to get all of this is the Apptainer/Singularity container
built from `../atac2grn.def` at the repo root:

```bash
# Build locally (from the repo root):
apptainer build atac2grn.sif atac2grn.def

# ...or pull the prebuilt image:
# https://drive.google.com/file/d/1UCvvgN3Bs85otx_AXmIazv1GznzPgU4C/view?usp=sharing

apptainer shell atac2grn.sif
```

### 2. Edit `config.yaml`

Inside the container (or on a workstation with the same tools on `PATH`),
copy and edit `Snakemake/config.yaml`:

```bash
cd Pipes/Snakemake
cp config.yaml my_run.yaml
```

Key fields (see inline comments in `config.yaml` for the full list):

| Key              | Meaning |
|------------------|---------|
| `input_dir`      | Folder containing one subfolder per sample; each subfolder is globbed for FASTQ files. |
| `output_dir`     | Where all pipeline output is written. |
| `mode`           | `"PE"` or `"SE"`. |
| `genome_build`   | Free-form label, also passed as `--organism` to `rgt-hint` (see caveat below). |
| `bowtie2_index`  | Prefix passed to `bowtie2 -x` (a prebuilt index, not a FASTA). |
| `peak_size`, `peak_min_dist` | HOMER `findPeaks` parameters. Default `peak_size: 50` matches the published hg19/mm9 pipeline that was validated against ChIP-seq; the chm13 pipeline used `500` but was never validated the same way, so treat `500` as a data point rather than a recommendation. |
| `mark_duplicates`| `true` runs Picard merge/sort/dedup before peak calling (matches the chm13 pipeline); `false` skips deduplication entirely (matches the original hg19/mm9 pipeline). |
| `fimo_bed`, `tss_bed` | Optional. Leave blank to skip the association-matrix (`bedops bedmap`) steps entirely. |

### 3. Run it

```bash
snakemake -s Snakefile --configfile my_run.yaml -n         # dry run: inspect the plan first
snakemake -s Snakefile --configfile my_run.yaml --cores 8  # real run
```

### Smoke test (no real data required)

`Snakemake/smoke_test/` contains a tiny synthetic paired-end fixture (a 4 kb
random reference + 40 matching read pairs) that lets you validate the whole
rule graph runs without any real reference genome or sequencing data. See
`Snakemake/smoke_test/README.md` for exact commands and its one caveat
(the footprinting step needs a registered "organism", same issue as chm13
below).

The Snakefile's rule graph (DAG structure, wildcard resolution, output
wiring) was also validated directly in this repo using
`snakemake -n` (Snakemake 9.27.0, local Python 3.13 install, no container)
against the smoke-test config in four combinations: PE+dedup, PE+no-dedup,
SE, and with `fimo_bed`/`tss_bed` set. All four produced the expected job
list with no errors. This validates pipeline wiring, not tool behavior:
actually executing the shell commands (bowtie2, Picard, HOMER, rgt-hint)
requires the container and was not performed in this environment (no
real/synthetic reference genomes were available to run bowtie2 itself
against; the "toy_genome.fa" fixture is a placeholder for exactly this).

## Using a genome RGT doesn't know about (e.g. chm13)

`rgt-hint` only ships built-in genome data for a fixed list of organisms:
hg19, hg38, mm9, mm10, mm39, tair10, zv9, zv10 (see the `data/` directory of
[CostaLab/reg-gen](https://github.com/CostaLab/reg-gen)). **chm13 is not on
this list.** If `genome_build` in your config is chm13 (or any other
unsupported build), `rgt-hint footprinting` will fail at that step unless you
first register it as a custom RGT organism:

1. RGT looks for organism definitions in a file at
   `$RGTDATA/data.config.user` (`$RGTDATA` defaults to `/root/rgtdata` inside
   the container, see `atac2grn.def`'s `%environment` block, or
   `$HOME/rgtdata` otherwise).
2. This repo's `Snakemake/data.config.user` already contains a `[chm13]`
   section (genome FASTA, chromosome sizes, gene models, etc.) from the
   original authors -- but Snakemake does not use it automatically. You must
   copy or merge it into RGT's data directory yourself, e.g.:
   ```bash
   mkdir -p "$RGTDATA"
   cp Pipes/Snakemake/data.config.user "$RGTDATA/data.config.user"
   ```
   (or merge its `[chm13]` section into an existing `data.config.user` if you
   already have one), and make sure the paths referenced inside that file
   (chm13 FASTA, chrom sizes, annotation files, etc.) actually exist on disk
   where you're running the pipeline -- they are not shipped with this repo
   or the container.
3. Bias-correction tables (Tn5/DNase cutting bias) are k-mer based and
   organism-agnostic, so no extra per-organism bias-table file is needed
   beyond the genome/annotation files above.

If you only need hg19, hg38, mm9, or mm10, none of this is necessary --
`genome_build` in `config.yaml` just needs to match one of those names.

## `fimo_bed`/`tss_bed` and genome coordinates

The original chm13 (`Snakefile-v4`) pipeline set `fimo_bed`/`tss_bed`
equivalents to **hg19** annotation files, even though reads were aligned to
**chm13**. If you set both of these while also using a non-hg19
`genome_build`, you are implicitly assuming a coordinate liftover was already
applied to your footprints before the `perform_fmap`/`perform_pmap` rules
run -- this pipeline does not perform that liftover for you. Leave
`fimo_bed`/`tss_bed` blank unless you have deliberately arranged for matching
coordinate systems.

## Legacy Snakefiles

The five Snakefiles this pipeline used to consist of
(`Snakefile-v2-hg19PE/SE`, `Snakefile-v2-mm9PE/SE`, `Snakefile-v4`) are
preserved unmodified under `Snakemake/legacy/` for exact reproduction of the
original paper's NIH HPC cluster runs. They are not runnable outside that
cluster (they use Lmod `module load` commands and Slurm `--cluster`
scheduling). See `Snakemake/legacy/README.md` for details. All new work
should use `Snakemake/Snakefile` + `Snakemake/config.yaml` instead.

