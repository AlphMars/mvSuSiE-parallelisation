# create the prior matrix for the SNP z-scores
arg = commandArgs(T)
zscores_path=arg[1]

suppressPackageStartupMessages({
    library(data.table)
    library(mvsusieR)
    })

zscores <- fread(zscores_path)
zscores_matr <- as.matrix(zscores[,-1])
rownames(zscores_matr) <- zscores$SNP
# remove rows with NA or Inf values 
zscores_matr <- zscores_matr[is.finite(rowSums(zscores_matr)),]

data   = mash_set_data(zscores_matr, Shat=1)
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)
U.pca = cov_pca(data,5,subset=strong)
U.ed = cov_ed(data, U.pca, subset=strong)

saveRDS(U.ed, file="prior_matrix_zscores.rds")

