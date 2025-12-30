#!/bin/bash
#SBATCH --partition=fat_rome
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --job-name=mvsusie_loop_with_prior
#SBATCH --output=/path/to/mvsusie/results/mvsusie_loop_with_prior.%A.out 
#SBATCH --error=/path/to/mvsusie/results/mvsusie_loop_with_prior.%A.err 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=XXX

lava_results_df="/path/to/all_pheno_pairs.all_loci.results.bivar.lava.fdr_sign.tsv" 
zscores="/path/to/zscore/matrix/zscores_matrix.tsv"
ldsc_output="/path/to/GenomicSEM/LDSCoutput.RData"
snp_info="/path/to/ref_genotype/g1000_eur_maf001.bim"
prior_zscores="path/to/zscore/matrix/prior_matrix_zscores.rds"
path/to/zscore/matrix/zscores_matrix.tsv

module load 2022
module load R/4.2.1-foss-2022a
    
Rscript mvsusie_loop_with_prior.R \
    "$lava_results_df" \
    "$zscores" \
    "$ldsc_output" \
    "$snp_info" \
    "$prior_zscores"
