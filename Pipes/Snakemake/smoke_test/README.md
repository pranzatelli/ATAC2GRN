# Smoke test fixture

This folder contains a tiny, fully synthetic paired-end dataset used to
validate that the consolidated `../Snakefile` executes its rule graph
correctly end-to-end, inside the container built from `atac2grn.def`. It is
**not** biological data and produces no scientifically meaningful footprints
or peaks -- it exists purely to catch pipeline-wiring bugs (missing rules,
wrong wildcards, broken shell quoting, dead-code outputs) without needing a
real reference genome or real sequencing data.

Contents:
- `generate_toy_data.py` -- regenerates everything below deterministically
  (fixed random seed). You do not need to run this; the outputs are already
  committed.
- `toy_genome.fa` -- a 4 kb random synthetic "chromosome" (`toy_chr1`).
- `data/sample1/sample1_R1.fastq.gz`, `sample1_R2.fastq.gz` -- 40 paired-end
  reads (75 bp), each an exact substring of `toy_genome.fa`, so a real
  bowtie2 alignment will succeed deterministically.
- `config.smoke.yaml` -- a ready-to-use Snakemake config pointing at the
  fixture above.

## Running the smoke test

From inside the container (see `../../README.md` "Quickstart"):

```bash
cd Pipes/Snakemake
bowtie2-build smoke_test/toy_genome.fa smoke_test/toy_genome
snakemake -s Snakefile --configfile smoke_test/config.smoke.yaml -n      # dry run
snakemake -s Snakefile --configfile smoke_test/config.smoke.yaml --cores 2
```

The dry run (`-n`) was also verified without a container, using a local
`pip install snakemake` under Python 3.13, for all four config combinations
(PE/dedup, PE/no-dedup, SE, and with `fimo_bed`/`tss_bed` set) -- see the
project's commit history / PR description for the exact `snakemake -n`
output. That check only validates the Snakemake rule graph (DAG structure,
wildcard resolution, output wiring); it does not execute bowtie2/Picard/
HOMER/RGT, since those tools are not installed outside the container.

## Caveat: the footprinting step needs a registered "organism"

`config.smoke.yaml` sets `genome_build: "toy"`, which (like chm13 -- see the
main README) is **not** a built-in RGT organism, so `rgt-hint footprinting`
will fail at that step unless you also register a minimal custom "toy"
organism with RGT (same mechanism as `data.config.user`, just pointing at
`toy_genome.fa` instead of a real genome). This is expected and is not a bug
in the Snakefile: the smoke test is primarily meant to validate the
alignment -> dedup -> HOMER peak-calling portion of the pipeline, which is
the majority of the DAG and does not require any organism registration.
If you want to exercise the footprinting rule too, either register a "toy"
organism, or point `bowtie2_index`/`genome_build` at a small real genome
that RGT already knows about (e.g. a subset of hg19) and adjust the fixture
accordingly.
