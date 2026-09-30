#!/bin/sh

module load snakemake
snakemake -s Snakefile-v4 --unlock --rerun-incomplete --cores 30
snakemake -s Snakefile-v4 -j100 -p -r --cores 30 --cluster "sbatch --time={cluster.time} --partition={cluster.partition} --mem={cluster.mem} --cpus-per-task={cluster.cpus}" --cluster-config cluster-v3.json --keep-going --latency-wait 120 --rerun-incomplete

