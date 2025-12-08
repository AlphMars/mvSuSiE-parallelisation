# --- mvsusie_loop.R --- #
# Run mvsusie_rss for each group of traits in a given locus.

#################################
# --- 1. ENVIRONMENT SET UP --- #
#################################
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
lava_results_path <- args[1]
zscores_matrix_path <- args[2]
ldsc_output_path <- args[3] 
snp_info_path <- args[4]
zscores_prior_matrix_path <- args[5]

# load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(parallel)
  library(mvsusieR)
  library(R.utils)
})

###############################
# --- 2. LOAD GLOBAL DATA --- #
###############################
cat("Loading global data...\n")

# Load LAVA results
lava_results_df <- fread(lava_results_path)

# Load LDSC results
load(ldsc_output_path)
cov_matrix_prior <- LDSCoutput_all_pheno_rev_ISI$S
residual_matrix <- cov2cor(LDSCoutput_all_pheno_rev_ISI$I)
rownames(residual_matrix) <- colnames(residual_matrix) <- rownames(cov_matrix_prior) <- colnames(cov_matrix_prior)

# Define trait information using named vectors
sample_sizes <- c(...)
  # specify the sample size of each trait in the exact same order of the traits below
  # when sample size varies across SNPs, calculate the median of the SNP-specific sample sizes
traits <- c("phen1", "phen2", ...) 
  # specify the exact names of the traits as they appear in the Z-scores matrix
  # and in the same order as the sample sizes above
names(sample_sizes) <- traits
rm(N, traits)

# create a variable with all the loci and, for each one of them, the associated list of traits
    locus_job_list <- lava_results_df %>%
    # Reshape the data from wide to long
    pivot_longer(
        cols = c(phen1, phen2),
        names_to = "phen_num",
        values_to = "trait"
    ) %>%
    # Group by each locus
    group_by(locus, chr, start, stop) %>%
    # Create a list of all unique traits for that locus
    summarise(
        traits_in_locus = list(unique(trait)),
        .groups = 'drop'
    )


################################################
# --- 3. DEFINE THE CORE ANALYSIS FUNCTION --- #
################################################
analyze_locus_pair <- function(i) {
  
  # --- Setup ---
  # Get task info from the i-th row of the global loci list
  # Note: assumes locus_job_list's column names to be: "locus", "chr", "start", "stop", "traits_in_locus"
  current_task <- locus_job_list[i, ]
  traits_to_analyze <- unlist(current_task$traits_in_locus)
  num_traits <- length(traits_to_analyze)
  locus_id <- current_task$locus
  locus_chr <- current_task$chr
  locus_start <- current_task$start
  locus_stop <- current_task$stop
  
  # Load SNP info
  snp_info <- fread(snp_info_path, col.names = c("chr", "snp_id", "cm", "pos", "a1", "a2"))

  # Define file names and core input files
  snp_list_file <- paste0("temp_snps_job_", i, ".txt")
  ld_out_file <- paste0("temp_ld_job_", i)  

  # Print a message when the job STARTS to track progress
  cat(paste0("Starting job ", i, ": Locus ", locus_id, " with ", num_traits, " traits.\n"))

  result <- tryCatch({

    #############################  
    # --- 4. PREPARE INPUTS --- #
    #############################
    
    # Step 1: Get all possible SNPs in the locus.
    locus_snps_initial <- snp_info %>%
      filter(chr == locus_chr, pos >= locus_start, pos <= locus_stop) %>%
      pull(snp_id)
    
    # Step 2: load the Z-matrix
    cols_to_select <- c("SNP", traits_to_analyze)
    locus_z_df <- fread(zscores_matrix_path, select = cols_to_select)
        
    # Step 3: Intersect with SNPs available in the Z-matrix and filter
    snps_in_locus <- intersect(locus_snps_initial, locus_z_df$SNP) 
    locus_z_df_subset <- locus_z_df[SNP %in% snps_in_locus, ]
    complete_cases_indices <- complete.cases(locus_z_df_subset[, ..traits_to_analyze])
    locus_z_df_filtered <- locus_z_df_subset[complete_cases_indices, ]

    if (nrow(locus_z_df_filtered) < 2) {stop("Skipping: Fewer than 2 common, non-NA SNPs.")}

    locus_z_matrix <- as.matrix(locus_z_df_filtered[, ..traits_to_analyze])
    final_snps_for_analysis <- locus_z_df_filtered$SNP # SNP IDs
    rownames(locus_z_matrix) <- final_snps_for_analysis # Set row names to SNP IDs

    if (length(final_snps_for_analysis) < 2) {stop("Skipping: Fewer than 2 common, non-NA SNPs in this locus.")}



    # Step 4: Use this clean list to generate the LD matrix.
    writeLines(final_snps_for_analysis, con = snp_list_file)
      
    plink_command <- paste(
      "~/plink1/plink",
      "--bfile /path/to/ref_genotype/g1000_eur_maf001",
      "--extract", snp_list_file,
      "--r square gz",
      "--threads 1",
      "--out", ld_out_file
    )
    # set "--threads 1" to control the number of threads used by PLINK and avoid conflicts with parallel execution
    system(plink_command, ignore.stdout = TRUE, ignore.stderr = TRUE)
      
    # Step 5: Process the LD matrix
    ld_matrix_raw <- fread(paste0(ld_out_file, ".ld.gz"), header = FALSE)
    ld_matrix <- as.matrix(ld_matrix_raw)
      
    # The order is guaranteed to match final_snps_for_analysis
    rownames(ld_matrix) <- final_snps_for_analysis
    colnames(ld_matrix) <- final_snps_for_analysis
      


    # Step 6: Create the prior matrices and the weights
    # For a multivariate analysis, we must provide a prior variance MATRIX.
    prior_matrices <- list()
    # covariance PCA matrices calculated from the SNPs zscores
    zscores_prior <- readRDS(zscores_prior_matrix_path)
    # identity matrix
    identity_prior <- diag(1, length(sample_sizes))
    rownames(identity_prior) <- colnames(identity_prior) <- colnames(cov_matrix_prior)
    # variance diagonal matrix
    cov_diag_prior <- matrix(0, length(sample_sizes), length(sample_sizes))
    diag(cov_diag_prior) <- diag(cov_matrix_prior)
    colnames(cov_diag_prior) <- rownames(cov_diag_prior) <- colnames(cov_matrix_prior)

    prior_matrices$covariance <- cov_matrix_prior[traits_to_analyze, traits_to_analyze]
    prior_matrices$identity <- identity_prior[traits_to_analyze, traits_to_analyze]
    prior_matrices$cov_diag <- cov_diag_prior[traits_to_analyze, traits_to_analyze]
    prior_matrices$zscores_ED_PCA_1 <- zscores_prior$ED_PCA_1[traits_to_analyze, traits_to_analyze]
    prior_matrices$zscores_ED_PCA_2 <- zscores_prior$ED_PCA_2[traits_to_analyze, traits_to_analyze]
    prior_matrices$zscores_ED_PCA_3 <- zscores_prior$ED_PCA_3[traits_to_analyze, traits_to_analyze]
    prior_matrices$zscores_ED_PCA_4 <- zscores_prior$ED_PCA_4[traits_to_analyze, traits_to_analyze]
    prior_matrices$zscores_ED_PCA_5 <- zscores_prior$ED_PCA_5[traits_to_analyze, traits_to_analyze]
    prior_matrices$zscores_ED_tPCA <- zscores_prior$ED_tPCA[traits_to_analyze, traits_to_analyze]

    # set fixed weights for the prior matrices
    num_matrix <- length(prior_matrices)
    mixture_weights <- rep(1/num_matrix, num_matrix)

    # Step 7: Create the prior object
    prior_mixture_object <- create_mixture_prior(
      list(
        matrices=prior_matrices,
        weights=mixture_weights
        ),
      null_weight=0
      )


    # Step 8: 
    # subset the residual variance matrix
    locus_residual_V <- residual_matrix[traits_to_analyze, traits_to_analyze]
    # Get the sample sizes for the traits
    locus_n <- sample_sizes[traits_to_analyze]
    

    ##############################
    # --- 5. RUN MVSUSIE_RSS --- #
    ##############################
    # --- Use a timeout to prevent long-running jobs from hanging indefinitely (NB: it's in seconds) --- #
    # NB: when `estimate_prior_variance = FALSE` each run should last less than 10 minutes.

    mvsusie_fit <- withTimeout({
      mvsusieR::mvsusie_rss(
        Z = locus_z_matrix,
        R = ld_matrix,
        N = min(locus_n), # Use the minimum sample size for the pair
        L = 10, # Try first with only 5 casual variants
        prior_variance = prior_mixture_object,
        residual_variance = locus_residual_V,
        estimate_prior_variance = FALSE,
        max_iter=1000,
        tol = 0.001
      )    
    }, timeout = 1200, onTimeout = "error")

      
      # --- RETURN A STRUCTURED RESULT ---
      list(
        status = "success",
        locus_id = locus_id,
        traits = traits_to_analyze,
        fit = mvsusie_fit
      )
      
    }, error = function(e) {
      list(
        status = "error",
        locus_id = locus_id,
        traits = traits_to_analyze,
        error_message = as.character(e)
      )
    }, finally = {
      # cleanup
      files_to_remove <- c(
        snp_list_file,
        paste0(ld_out_file, ".ld.gz"),
        paste0(ld_out_file, ".log"),
        paste0(ld_out_file, ".nosex")
      )
      for (f in files_to_remove) {
        if (file.exists(f)) {
          file.remove(f)
        }
      }

  })
    
  return(result)
}


##################################
# --- 6. EXECUTE IN PARALLEL --- #
##################################

# --- Get number of cores from the SLURM environment variable --- #
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
if (is.na(num_cores)) {num_cores <- 1}

cat("Starting analysis for", nrow(lava_results_df), "locus-pair combinations.\n")

# --- Create a socket cluster --- #
cl <- makeCluster(num_cores)
cat("Created a cluster with", length(cl), "workers.\n")

# --- Export the necessary global objects and functions to each worker node. --- #
clusterExport(cl, varlist = c(
  "locus_job_list",
  "snp_info_path",
  "zscores_matrix_path",
  "zscores_prior_matrix_path",
  "sample_sizes",
  "cov_matrix_prior",
  "residual_matrix",
  "analyze_locus_pair"
))

# --- Load the required packages on each worker node. --- #
clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(tidyr)
    library(mvsusieR)
    library(R.utils)
  })
})

cat("Exported data and loaded packages on all workers. Starting parLapply...\n")

# --- Apply the function in parallel using parLapply (safer alternative to mclapply). --- #
final_results_list <- parLapply(cl, 1:nrow(locus_job_list), analyze_locus_pair)

# --- Stop the cluster to free up resources. --- #
stopCluster(cl)
cat("Analysis complete! Cluster stopped. Saving results...\n")

# --- Save the final list of results --- #
saveRDS(final_results_list, file = "mvsusie_results_with_prior.rds")

cat("Done.\n")


# --- END OF SCRIPT --- #