#!/usr/bin/env python3
"""Generate a tiny, fully synthetic paired-end ATAC-seq-like smoke-test fixture.

This creates a small reference genome (toy_genome.fa) and paired-end reads
(data/sample1/sample1_R1.fastq.gz / sample1_R2.fastq.gz) that are exact
substrings of the reference, so a real bowtie2 alignment against a bowtie2
index built from toy_genome.fa will succeed deterministically. Read counts and
positions are deliberately small/toy: this is meant to validate that the
consolidated Snakefile's rules execute end-to-end inside the container (see
Pipes/README.md "Smoke test"), not to produce biologically meaningful
footprints or peaks.

Re-running this script regenerates the exact same output (fixed random seed),
so the generated files are also committed to the repo for convenience.
"""
import gzip
import random
from pathlib import Path

SEED = 42
GENOME_LEN = 4000
READ_LEN = 75
FRAGMENT_LEN = 200
N_PAIRS = 40

HERE = Path(__file__).resolve().parent


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def main() -> None:
    rng = random.Random(SEED)
    genome = "".join(rng.choice("ACGT") for _ in range(GENOME_LEN))

    fasta_path = HERE / "toy_genome.fa"
    with open(fasta_path, "w", newline="\n") as fh:
        fh.write(">toy_chr1\n")
        for i in range(0, len(genome), 70):
            fh.write(genome[i : i + 70] + "\n")

    sample_dir = HERE / "data" / "sample1"
    sample_dir.mkdir(parents=True, exist_ok=True)
    r1_path = sample_dir / "sample1_R1.fastq.gz"
    r2_path = sample_dir / "sample1_R2.fastq.gz"

    max_start = GENOME_LEN - FRAGMENT_LEN
    with gzip.open(r1_path, "wt", newline="\n") as r1_fh, gzip.open(
        r2_path, "wt", newline="\n"
    ) as r2_fh:
        for i in range(N_PAIRS):
            start = rng.randint(0, max_start)
            fragment = genome[start : start + FRAGMENT_LEN]
            read1 = fragment[:READ_LEN]
            read2 = revcomp(fragment[-READ_LEN:])
            qual = "I" * READ_LEN
            name = f"@toy_pair_{i}"
            r1_fh.write(f"{name}/1\n{read1}\n+\n{qual}\n")
            r2_fh.write(f"{name}/2\n{read2}\n+\n{qual}\n")

    print(f"Wrote {fasta_path}")
    print(f"Wrote {r1_path} and {r2_path} ({N_PAIRS} pairs)")


if __name__ == "__main__":
    main()
