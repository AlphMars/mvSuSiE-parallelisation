#!/bin/bash
#SBATCH --partition=rome
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --job-name=prior_zscores

zscores="/path/to/zscore/matrix/zscores_matrix.tsv"

module load 2022
module load R/4.2.1-foss-2022a

# zscores
Rscript prior_matrix_zscores.R "$zscores"
