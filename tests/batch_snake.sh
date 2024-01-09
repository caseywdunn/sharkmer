#!/bin/bash

#SBATCH --job-name=snake
#SBATCH --output=snake.txt
#SBATCH --time=1-00:00:00
#SBATCH --partition=ycga
#SBATCH --nodes=1                    # number of cores and nodes
#SBATCH --cpus-per-task=8           # number of cores
#SBATCH --mem-per-cpu=8G            # shared memory, scaling with CPU request

# For sharkmer, use --partition=ycga_bigmem    --cpus-per-task=16    --mem-per-cpu=60G
# For most other steups, use --partition=pi_dunn   --cpus-per-task=16   --mem-per-cpu=4G

module load SRA-Toolkit

module load miniconda
# conda create -n snakemake -c conda-forge -c bioconda genomescope2 jellyfish snakemake matplotlib pandas
conda activate snakemake

snakemake --cores $SLURM_CPUS_PER_TASK