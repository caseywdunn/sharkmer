#!/bin/bash

#SBATCH --job-name=snake
#SBATCH --output=snake.txt
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --nodes=1                    # number of cores and nodes
#SBATCH --cpus-per-task=8           # number of cores
#SBATCH --mem-per-cpu=18G            # shared memory, scaling with CPU request

module load miniconda
## conda create -n snakemake -c conda-forge -c bioconda genomescope2 jellyfish snakemake matplotlib pandas sra-tools
ml FastQC
conda activate snakemake

snakemake --touch --rerun-incomplete
snakemake --cores $SLURM_CPUS_PER_TASK --rerun-incomplete
