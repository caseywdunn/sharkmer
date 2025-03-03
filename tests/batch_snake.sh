#!/bin/bash

#SBATCH --job-name=snake
#SBATCH --output=snake.txt
#SBATCH --time=1-00:00:00
#SBATCH --partition=ycga
#SBATCH --nodes=1                    # number of cores and nodes
#SBATCH --cpus-per-task=24           # number of cores
#SBATCH --mem-per-cpu=18G            # shared memory, scaling with CPU request

# For sharkmer, use --partition=ycga_bigmem    --cpus-per-task=16    --mem-per-cpu=60G
# For most other steups, use --partition=pi_dunn   --cpus-per-task=16   --mem-per-cpu=4G

module purge
#module load SRA-Toolkit # Python from this module is conflicting with miniconda
#ml snakemake
#snakemake --cores 8 --rerun-incomplete 

module load miniconda
## conda create -n snakemake -c conda-forge -c bioconda genomescope2 jellyfish snakemake matplotlib pandas
ml FastQC
conda activate snakemake

snakemake --touch --cores 1 --rerun-incomplete
snakemake --cores $SLURM_CPUS_PER_TASK --rerun-incomplete
