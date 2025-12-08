#!/bin/bash
#SBATCH --partition=fat_rome
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --job-name=mvsusie_loop_with_prior
#SBATCH --output=/home/amartone/colocal/mvsusie_loop_with_prior.%A.out 
#SBATCH --error=/home/amartone/colocal/mvsusie_loop_with_prior.%A.err 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alfonso.martone@radboudumc.nl

lava_results_df="/home/amartone/colocal/all_pheno_pairs.all_loci.results.bivar.lava.fdr_sign.tsv"
zscores="/home/amartone/colocal/zscores_matrix.tsv"
ldsc_output="/home/amartone/genomicsem/genomicsem_cog_IR_all/LDSCoutput_all_pheno_rev_ISI.RData"
snp_info="/home/amartone/lava/ref_genotype/g1000_eur_maf001.bim"
prior_zscores="/home/amartone/colocal/prior_matrix_zscores.rds"

module load 2022
module load R/4.2.1-foss-2022a
    
Rscript mvsusie_loop_with_prior.R \
    "$lava_results_df" \
    "$zscores" \
    "$ldsc_output" \
    "$snp_info" \
    "$prior_zscores"
