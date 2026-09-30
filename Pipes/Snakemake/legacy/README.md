# Legacy Snakefiles

The five files in this folder (`Snakefile-v2-hg19PE`, `Snakefile-v2-hg19SE`,
`Snakefile-v2-mm9PE`, `Snakefile-v2-mm9SE`, `Snakefile-v4`) plus
`cluster-v3.json` and `run_snake.sh` are kept **unmodified** exactly as they
were used to produce the results in the original paper, run on the authors'
NIH HPC cluster.

They are **not runnable as-is** outside that cluster:
- They use Lmod `module load bowtie|picard|homer|bedtools|rgt|bedops` commands
  that only exist on that cluster's environment modules system.
- `Snakefile-v4` and `run_snake.sh`/`cluster-v3.json` assume a Slurm scheduler
  (`sbatch`) is available via `snakemake --cluster`.
- Paths (`genomes/hg19`, `CHM13/Bowtie2Index/chm13`, `hg19/v3_fimo.bed`, etc.)
  are hardcoded to locations on that cluster.

They are kept here purely for provenance/exact reproduction of the published
analysis. **For actually running the pipeline** (in the Apptainer container
built from `atac2grn.def`, or on a plain workstation with the same tools on
PATH), use `../Snakefile` together with `../config.yaml` instead — see
`Pipes/README.md` for the Quickstart. That consolidated Snakefile documents,
rule by rule, which real logic differences between these five legacy files it
preserves (via config.yaml options) versus which it deliberately changes (and
why) compared to the legacy behavior.
