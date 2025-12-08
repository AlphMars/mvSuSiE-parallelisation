# mvSuSiE-parallelization
Script and README files for mvSuSiE parallelization across multiple traits and genomic chunks and to extract summary informations

The multivariate extension of the Sum of Single Effects (mvSuSiE) model performs fine-mapping on given genomic regions allowing to identify high-probability causal variants from GWAS summary statistics. It learns the patterns of shared genetic effects among traits and adapts automatically to the pattern. It's also possible to perform colocalization analyses. See the [paper][1] and the [GitHub][2] for more.
It's an extremely time- and resource-efficient tool, however for large number of traits parallelization is needed. I created a Slurm job for HPC clusters that can be used for this purpose, integrating available information from the [original tutorial][3] and personal communications with the authors.


### Inferences
- PIP (posterior inclusion probability)
- CSs (credibl sets)
  - designed to capture with high probability, at least one causal SNP
  - informally, each one represents an independent association signal in the data
  - size (number of SNPs) ~ how precisely one can pinpoint the causal SNP underlying this signal
- separate inferences in two questions
    1. which SNPs are causal for at least 1 trait
       - **cross-trait PIPs and CSs**
       - summarize the inferences across all traits
       - >threshold at >0.95
    2. which traits does each identified CS affect?
       - **lfsr** (local false sign rate)
       - trait-wise measure, for each SNP in each trait
       - >threshold at <0.01
       - **average lfsr** as a trait-wise measure (SNPs in the same CS are likely in high LD)

<br>

## Things you need to perform mvSuSiE
Being J the number of SNPs and T the number of traits:

### 1. Z-score matrix $Z$
A J x T matrix of Z-scores. You can calculate it from the munged sumstats:
```r
library(data.table)
files <- c(...) # list your files, eventually with full or relative path
traits <- c(...) # list the traits, keeping the same order of the files
mylist <- list() # create an empty list to store the extracted data

# loop through the files
# if you're using munged sumstats, they're likely to be in the '.sumstats.gz' format
# you can use fread() to read those files
for (i in 1:length(files)) {
    df <- fread(files[i]) %>% as.data.frame()   # convert in dataframe
    trait <- traits[i]
    clean_df <- df %>% select(SNP, Z)   # extract the needed columns
    colnames(clean_df)[2] <- trait      # Rename the Z-score column to be trait-specific
    mylist[[i]] <- clean_df    # Add the cleaned data frame to the list
    }

combined_df <- reduce(mylist, full_join, by = "SNP")    # Combine all data frames

# Set the SNP column as the row names and convert to a matrix
rownames(combined_df) <- combined_df$SNP
z_matrix <- combined_df %>% 
select(-SNP) %>%
as.matrix()
write.table(z_matrix, "zscores_matrix.tsv", col.names=T, row.names=T, sep="\t", quote=F)
```

### 2. Sample size vector $n$
A vector of sample sizes, with length T, in the same order as the columns of the Z-score matrix.
It must be added in the R script

### 3. LD matrix $R$
A J x J LD matrix. This will be calculated directly in the script using PLINK and a suitable reference panel (e.g., 1000 Genomes).
The rows and columns of R must be in the exact same order as the rows (SNPs) of $Z$ 

### 4. Prior_variance
Controls the expected size of causal effects. This is the trickiest part: you need to provide several matrices representing possible alternative conditions. Part of them are estimated from the z-score matrix $Z$ using the `prior_matrix_zscores.sh` script. Another part will be extracted from the LDSC results directly in the final script.
I came to this conclusion after contacting the [authors][4], they say that using only the LDSC matrices as prior wouldn't be sufficient since it's too strict. mvSuSiE works better if a set of different matrices is provided with different hypothesis, like:
- identity matrix (diagonal=1, off diagonal = 0): each trait has indipendent effects
- single-trait variance matrices
- zscore matrices

### 5. Residual_variance
It estimates the prior belief about the noise. The intercept matrix obtained from LDSC is a valid estimate.

### Tips and notes
In case you ran the LDSC through the GenomicSEM package, you will find there:
- covariance matrices in `LDSCoutput.RData$S`
- correlation matrices in `LDSCoutput.RData$S_Stand`
- sample size matrix in `LDSCoutput.RData$N`

>Of course, before starting you need to install the mvSuSiE package by following the instructon in the [GitHub][2] repository.

<br>

Running mvSuSiE
--------------------------------------------

### 1. Estimate z-score prior matrices
This step is done by calling the `prior_matrix_zscores.sh` script to which you should provide the path to your z-score matrix. The job calls the `prior_matrix_zscores.R` script which performs a Multivariate Adaptive Shrinkage (mash) analysis, extimating the covariance structures of the z-score matrix using Principal Component Analysis and Extreme Deconvolution.

### 2. Run mvSuSiE
To run mvSuSiE you must first adapt the two scripts. 
In `mvsusie_loop.sh`, you must specify:
- the lava results dataframe or more generally a **list of chromosome regions per trait** that will be used to run mvSuSiE on each region separatedly. The structure of the dataframe must be like the following:
```
    locus   chr     start   end     phenotype
    9       1       6136815 7711794 B_fluid_int
    17      1       15755231        16732168        B_fluid_int
    58      1       67761891        68633860        B_fluid_int
    62      1       72513120        73992170        B_fluid_int
    70      1       80781443        82123368        B_fluid_int
    155     1       201067953       202583884       B_fluid_int
    230     2       26894103        28819510        B_fluid_int
    249     2       44104111        45189468        B_fluid_int
    266     2       57952946        59251996        B_fluid_int
```
- the **z-score matrix**
- the **LDSC output**, an object including the covariance matrix and the intercept matrix
- the **SNP reference panel**, e.g. *g1000_eur_maf001.bim* from the 1000 genome project
- the **prior z-score matrices** calculated above

This job will cal the `mvsusie_loop.R` script, in which you must further specify:
- the **sample size vector**
- the **trait vector**
- the **path to reference genotype** for the `--bfile` inside the plink command (around line 128)

<br>

Checking the results
------------------------------------
Results from mvSuSiE are dense of informations. Here's what I checked for the purpose of my paper.

### Successful and unsuccessful fits, null CS 
```r
res <- readRDS("mvsusie_results_with_prior.rds")
library(purrr)
library(dplyr)
library(data.table)


# 1. ---< filter to separate successful and unsuccessful fits >--- #
    unsuccessful_fits <- keep(res, ~ .$status != "success")
    successful_fits <- keep(res, ~ .$status == "success")

    # Check how many succeeded
    cat("Number of successful fits:", length(successful_fits), "\n")
        # Number of successful fits: 219

# 2. ---< extract data about the unsuccessful >--- #
    failed_trait_pairs <- map(unsuccessful_fits, ~ .$traits)
    failed_locus_ids <- map_int(unsuccessful_fits, ~ .$locus_id)
    error_messages <- map_chr(unsuccessful_fits, ~ .$error_message)
    # combine data into a useful data frame for review
    error_summary <- tibble(
      locus_id = failed_locus_ids,
      trait_pair = failed_trait_pairs,
      error = error_messages
    )

    print(error_summary)
        # # A tibble: 2 × 3
        #   locus_id trait_pair error                                                     
        #      <int> <list>     <chr>                                                     
        # 1      953 <chr [3]>  "Error in doTryCatch(return(expr), name, parentenv, handl…
        # 2      959 <chr [4]>  "Error in doTryCatch(return(expr), name, parentenv, handl…
        
    res_summary <- tibble(
        locus_id = map_int(res, "locus_id"),
        traits = map_chr(res, ~ paste(.$traits, collapse=", ")),
        pip = map(res, ~ .$fit$pip),
        mean_pip = map_dbl(res, ~ mean(.$fit$pip)),
        lfsr = map(res, ~ .$fit$lfsr),
        mean_lfsr = map_dbl(res, ~ mean(.$fit$lfsr)),
        cs = map(res, ~ .$fit$sets$cs)
    )


# 4. --- subset based on the CS --- #
    res_summary_null_cs <- filter(res_summary, map_lgl(cs, is.null))
    # # or
    # res_summary_null_cs <- res_summary %>%
    #     filter(map_lgl(cs, is.null))
      # 1. map_lgl(cs, is.null) goes through each item in the cs column.
      # 2. For the first row, it computes is.null(<dbl [286]>) which is FALSE.
      # 3. For the second row, it computes is.null(<
      # 4. If it found a row where cs was NULL, is.null(<NULL>) would be TRUE.
      # 5. It returns a logical vector like [TRUE, FALSE, FALSE, TRUE, ...].
      # 6. filter() then uses this logical vector to keep only the rows corresponding to TRUE.
    
    res_null_cs <- keep(res, ~is.null(.x$fit$sets$cs))
      # 1. keep() iterates through each top-level element of the res list.
      # 2. In each iteration, the current element is referred to as .x.
        # The ~ and .x are a special shorthand for creating an anonymous function
        # It's equivalent to: function(x) is.null(x$fit$sets$cs)
      # 3. The expression is.null(.x$fit$sets$cs) is evaluated.
      # 4. If it's TRUE, keep() saves that element (.x) for the final output list. If FALSE, it's discarded.
    res_cs <- discard(res, ~is.null(.x$fit$sets$cs))
```

<br>

### Get PIP, lfsr, and loci
```r
  traits <- map(res, ~ .$traits)
  traits_length <- map_int(res, ~ length(.$traits))
  locus_ids <- map_int(res, ~ .$locus_id)
  traits_per_locus <- tibble(
    locus_id=locus_ids, 
    number_of_traits=traits_length
    )
  
# ---< extract, for each cs in each locus, its pip and lfsr >--- #
  # get locus names
    locus_id <- map(res, ~ .$locus_id)
    loci <- unlist(locus_id)
    locus_name <- paste0("Locus ", loci)
  # create the summary list  
    traits <- map(res, ~ .$traits)
    # get cs list for each locus
      cs_vector <- list()
      for (i in 1:length(loci)) {
        if (!is.null(res[[i]]$fit$sets$cs)) {
          cs_vector[[i]] <- unlist(res[[i]]$fit$sets$cs)
        } else {next}
      }
      names(cs_vector) <- locus_name

    # get pip and lfsr for each CS
      pip <- list()
      lfsr <- list()
      for (i in 1:length(loci)) {
        if (!is.null(res[[i]]$fit$sets$cs)) {
          pip[[i]] <- res[[i]]$fit$pip[cs_vector[[i]]]
          lfsr[[i]] <- res[[i]]$fit$lfsr[cs_vector[[i]],]
          colnames(lfsr[[i]]) <- traits[[i]]
        } else {next}
      }     
      names(pip) <- locus_name
      names(lfsr) <- locus_name
      
      # N.B.: while there is only one pip for each SNP, you get one lfsr for each SNP and each trait, that's y additional steps are needed
    
    # combine everything in one tible for each locus, plus trait names and locus id
      cs_summary <- list()
      length(cs_summary) <- length(loci)
      names(cs_summary) <- locus_name
      for (i in 1:length(loci)) {
        if (!is.null(res[[i]]$fit$sets$cs)) {
          cs_summary[[i]]$locus_id <- locus_id[[i]]
          cs_summary[[i]]$traits <- traits[[i]]
          cs_summary[[i]]$cs <- tibble(
            cs = names(cs_vector[[i]]),
            index = cs_vector[[i]],
            SNP = names(pip[[i]]),
            pip = pip[[i]],
            !!!as_tibble(lfsr[[i]])
          )
        } else {
          cs_summary[[i]]$locus_id <- locus_id[[i]]
          cs_summary[[i]]$traits <- traits[[i]]
          }
      }

# ---< extract loci involved in at least 3 traits >--- #
    cs_summary_filtered <- keep(cs_summary, ~ length(.$traits) >= 3)
    
# ---< extract the CS from these >--- #
    var <- unlist(map(cs_summary_filtered, ~ .$cs$SNP)) %>% sort() %>% unique()
```

<br>

### Filter for PIP > 0.95
```r
# ---< Use map() to iterate through each element of cs_summary >--- #
    cs_pip_95 <- map(cs_summary, function(locus_result) {
    
        # Check if the locus has a 'cs' tibble to filter
        if (!is.null(locus_result$cs)) {
            # Use dplyr's filter to keep only the high-PIP rows
            locus_result$cs <- locus_result$cs %>%
            filter(pip >= 0.95)
            }
        
        # return the modified locus_result (with the filtered or original cs tibble)
        return(locus_result)
    })

# ---< remove loci that now have an empty cs tibble >--- #
    cs_pip_95_final <- keep(cs_pip_95, ~ !is.null(.$cs) && nrow(.$cs) > 0)
    saveRDS(cs_pip_95_final, "mvsusie_results_with_prior_CS.pip_95.rds")
    
# ---< extract loci with at least three involved traits >--- #
    cs_pip_95_final_3 <- keep(cs_pip_95_final, ~ length(.$traits) >= 3)
    saveRDS(cs_pip_95_final_3, "mvsusie_results_with_prior_CS.pip_95_traits_3.rds")

# ---< extract the CS from these and create the table for SNPnexus >--- #
    var_all <- unlist(map(cs_pip_95, ~ .$cs$SNP)) %>% sort() %>% unique()
    dbsnp_all <- paste0(rep("dbsnp", length(var_all)))
    variants_all <- tibble(a=dbsnp_all, snp=var_all)
    write.table(variants_all, "mvsusie_results_with_prior.variants.pip_95_alltraits.tsv", col.names = F, row.names = F, quote = F, sep=" ")

    var <- unlist(map(cs_pip_95_final_3, ~ .$cs$SNP)) %>% sort() %>% unique()
    dbsnp <- paste0(rep("dbsnp", length(var)))
    variants <- tibble(a=dbsnp, snp=var)
    write.table(variants, "mvsusie_results_with_prior.variants.pip_95_traits_3.tsv", col.names = F, row.names = F, quote = F, sep=" ")
```

<br>

### Gather results in a table and filter for lfsr < 0.01
```r
library(purrr)
library(dplyr)
library(tidyverse) # Good to have for data manipulation
res95 <- readRDS("mvsusie_results_with_prior_CS.pip_95.rds")
# in case you want to rename/reorder your traits you can use a named vector
traits <- c("phen1", "phen2", ...)  # define your traits in the desired order
load("/path/to/GenomicSEM/LDSCoutput.RData") 
colnames <- colnames(LDSCoutput$S)
names(colnames) <- traits   # the order of the traits should match, in case not, reorder it

# ---< Step 1: Process each locus to prepare it for binding >---
  # use map() to iterate through each element of your 'res95' list.
  # for each locus, add the new ID columns to its 'cs' tibble.
  ###### WARNING : this code assumes that res95 is a named list, e.g., names(res95) are "Locus 9", "Locus 57", etc.

  processed_list <- map(res95, function(locus_result) {
  
    # handle cases where a locus might have a high-PIP SNP but the cs table is empty/null
    if (is.null(locus_result$cs) || nrow(locus_result$cs) == 0) {
      return(NULL)
    }
  
  # get the locus ID from within the list element
  loc_id <- locus_result$locus_id
  
  # modify the 'cs' tibble to add the new ID column
  locus_result$cs %>%
    mutate(
      locus_cs_id = paste0("Locus_", loc_id, "_", cs),
      .before = 1 # A nice trick to put the new column at the front
    )
  })

# ---< Step 2: Combine the list of tibbles into a single data frame >---
  # bind_rows() aligns columns by name and fills missing values with NA.
  final_results_df <- bind_rows(processed_list)

# ---< Step 3: Select and reorder columns into your desired final format >---
  final_results_df <- final_results_df %>%
    select(
      locus_cs_id,
      index,
      SNP,
      pip,
      everything(),  # This gets all other columns (the traits)
      -cs,            # We can remove the original 'cs' column as it's now in locus_cs_id
    )

  final_results <- final_results_df %>% rename(!!!colnames)
  # the !!! (pronounced "bang-bang-bang") operator is from the rlang package (part of the Tidyverse). 
  # It tells rename() to "unquote and splice" the named vector.
  # This step is ment for renaming the traits, you can skip it by simply uncommenting the following command:
  # final_results <- final_results_df 
  
  final_results_reordered <- final_results %>%
  select(
    # First, list your fixed ID columns
    locus_cs_id,
    SNP,
    pip,
    # use any_of() to select the trait columns that exist in your
    # data frame, in the order specified by your template vector.
    any_of(traits)
  )

  write.table(final_results, "mvsusie_results_with_prior_CS.pip_95.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# ---< Step 4: Filter for lfsr <0.01 >---
  filtered_lfsr_df <- final_results %>%
    mutate(
      # The across() function applies a function to a selection of columns
      across(
        # 1. Select all columns EXCEPT the ones listed here
        !c(locus_cs_id, index, SNP, pip),
        
        # 2. Apply this function to each selected column
        ~ if_else(.x >= 0.01, NA_real_, .x)
      )
    )
  
  write.table(filtered_lfsr_df, "mvsusie_results_with_prior_CS.pip_95.lfsr_01.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# ---< Step 5 : keep rows with at least one non-NA value >---
  final_significant_results <- filtered_lfsr_df %>%
    filter(
      # if_any() returns TRUE for a row if the condition is met in any of the specified columns
      if_any(
        # 1. Specify the columns to check: all columns EXCEPT the ID columns
        !c(locus_cs_id, index, SNP, pip),
        
        # 2. Specify the condition: is the value NOT NA?
        # The ~ creates a mini-function, .x is the value in the cell
        ~ !is.na(.x)
      )
    )
  write.table(final_significant_results, "mvsusie_results_with_prior_CS.pip_95.lfsr_01.noNA.tsv", sep="\t", row.names = F, col.names = T, quote = F)
```

The `final_significant_results` variable should look like this:
```
# A tibble: 847 × 25
   locus_cs_id index SNP         pip    Fluid      WHR   BMI Pairs meanRT `TMT-A`  MetS   PAL   Num   FPG   SDS  T2DM   FPI Tower HbA1c
   <chr>       <int> <chr>     <dbl>    <dbl>    <dbl> <dbl> <dbl>  <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
 1 Locus_9_L1    337 rs200448  0.999 9.43e- 4 9.43e- 4    NA    NA     NA      NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
 2 Locus_9_L2    394 rs277679  1.000 5.20e- 5 5.20e- 5    NA    NA     NA      NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
 3 Locus_9_L3    334 rs200440  1     1   e-20 1   e-20    NA    NA     NA      NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
 4 Locus_9_L4    398 rs278020  1     1   e-20 1   e-20    NA    NA     NA      NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
 5 Locus_9_L5    418 rs350704… 1     1   e-20 1   e-20    NA    NA     NA      NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
# ℹ 842 more rows
# ℹ 6 more variables: IFC <dbl>, ISI <dbl>, `2hGlu` <dbl>, Matrix <dbl>, Obesity <dbl>, `TMT-B` <dbl>
# ℹ Use `print(n = ...)` to see more rows
```
It lists the PIP for each SNP and the corresponding lfsr values for that SNP with the corresponding traits.

<br>

### Merging with LAVA
If you ran LAVA or another kind of local correlation analysis or you whish to map back each variant to its genomic region, you can provide it and merge it with the results above.
```r
lava <- read.table("/path/to/lava_results/all_pheno_pairs.all_loci.results.bivar.lava.fdr_sign.renamed.tsv", sep="\t", header=T)
    # a table with the following structure:
        # locus   chr     start   end     phenotype
        # 9       1       6136815 7711794 B_fluid_int
        # 17      1       15755231        16732168        B_fluid_int
        # 58      1       67761891        68633860        B_fluid_int
        # 62      1       72513120        73992170        B_fluid_int
        # 70      1       80781443        82123368        B_fluid_int
        # 155     1       201067953       202583884       B_fluid_int
        # 230     2       26894103        28819510        B_fluid_int
        # 249     2       44104111        45189468        B_fluid_int
        # 266     2       57952946        59251996        B_fluid_int

# ---< Step 0: create locus names in the forman "chrN:start-stop" >---
    lava_renamed <- lava %>% mutate(locus_name=paste0("chr", lava$chr, ":", lava$start, "-", lava$stop), locus_id=paste0("Locus_", locus))
    lava_renamed <- lava_renamed %>% 
        select(locus_id, locus_name) %>%
        distinct(locus_id, .keep_all = TRUE)    # the original file may have duplicated locus (the hotspots)

# ---< Step 1: Process each locus to prepare it for binding >---
  # use map() to iterate through each element of your 'res95' list.
  # for each locus, add the new ID columns to its 'cs' tibble.

  processed_list <- map(res95, function(locus_result) {
  
    # handle cases where a locus might have a high-PIP SNP but the cs table is empty/null
    if (is.null(locus_result$cs) || nrow(locus_result$cs) == 0) {return(NULL)}
  
    loc_id <- locus_result$locus_id    # get the locus ID from within the list element
  
    # modify the 'cs' tibble to add the new ID column
    locus_result$cs %>%
        mutate(
        locus_cs_id = paste0("Locus_", loc_id, "_", cs),
        .before = 1 # put the new column at the front
        )
  })

# ---< Step 2: Combine the list of tibbles into a single data frame, add the locus informations, rename and reorder >---

  final_results_df <- bind_rows(processed_list)

  final_results_loci <- final_results_df %>% 
    separate(locus_cs_id, into = c("locus_prefix", "locus_id", "cs_id"), sep = "_", remove = FALSE) %>%
    mutate(locus_id = paste0(locus_prefix, "_", locus_id)) %>% 
    rename(!!!colnames) %>%
    # grab the column you need
    select(
        locus_cs_id,
        locus_id,
        index,
        SNP,
        pip,
        any_of(traits),  # use any_of() to select the trait columns in the order specified by your template vector.
        -cs
        )


# ---< Step 3: Add the locus coordinates >---
  final_results_loci_coord <- final_results_loci %>% left_join(lava_renamed, by = "locus_id") %>%
      select(
        locus_name, # The new position column
        locus_cs_id,
        index,
        SNP,
        pip,
        everything(),   # All other columns (traits, etc.)
        -locus_id
    )

  write.table(final_results_loci_coord, "mvsusie_results_with_prior_CS.pip_95.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# ---< Step 4: Filter for lfsr<0.01 >---
  filtered_lfsr_df <- final_results_loci_coord %>%
    mutate(
      # The across() function applies a function to a selection of columns
      across(
        # select all columns EXCEPT the ones listed here
        !c(locus_name, locus_cs_id, index, SNP, pip),
        # apply the function to each selected column
        ~ if_else(.x >= 0.01, NA_real_, .x)
      )
    )
  
  write.table(filtered_lfsr_df, "mvsusie_results_with_prior_CS.pip_95.lfsr_01.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# ---< Step 5 : keep rows with at least one non-NA value >---
  final_significant_results <- filtered_lfsr_df %>%
    filter(
      if_any(
        !c(locus_name, locus_cs_id, index, SNP, pip),
        ~ !is.na(.x)
      )
    )
  write.table(final_significant_results, "mvsusie_results_with_prior_CS.pip_95.lfsr_01.noNA.tsv", sep="\t", row.names = F, col.names = T, quote = F)
```

<br>

### Extract general informations
According to how results were reported in Zou 2023.

```r
CS <- cs_summary

# ---< Step 1 : how many CS per candidate locus >---
    CS_count <-  map(CS,~ nrow(.$cs)) %>% bind_rows() %>% t()
    summary(CS_count)
        #        V1        
        #  Min.   : 2.000  
        #  1st Qu.: 3.000  
        #  Median : 6.000  
        #  Mean   : 5.756  
        #  3rd Qu.: 8.000  
        #  Max.   :12.000  

# ---< Step 2 : how many SNPs per CS >---
    res <- readRDS("mvsusie_results_with_prior.rds")

# 1. filter for successful runs.
    successful_fits <- keep(res, ~ .$status == "success" && !is.null(.$fit$sets$cs))
    null_cs <- keep(res, ~ is.null(.$fit$sets$cs))
    errors <- keep(res, ~ .$status == "error")

# 2. summarize the umber of CS per locus
    # The map_dfr() function is a shortcut for running map()
    # and then bind_rows() on the results.
    cs_size_summary <- map_dfr(successful_fits, function(locus_result) {
        loc_id <- locus_result$locus_id    # get the locus ID for this result
        cs_list <- locus_result$fit$sets$cs    # extract the nested list of credible sets
        if (length(cs_list) == 0) {return(NULL)}    # handle the rare case where cs exists but is empty
        cs_lengths <- map_int(cs_list, length)    # get the length of each credible set
        tibble(
            cs_id = names(cs_lengths),
            num_variants = cs_lengths,
            locus_id = loc_id # Add the locus ID to every row
        )
    }) %>%  
        # reorder columns and create the combined ID
        mutate(locus_cs_id = paste0("Locus_", locus_id, "_", cs_id), .before = 1) %>%
        mutate(locus_seq = paste0("Locus_", locus_id), .before =1) %>%
        select(-locus_id, -cs_id) 

# 3. get the toatl number of CS 
    nrow(cs_size_summary)

# 4. get the number of CS for each locus by counting duplicates
    cs_count <- cs_size_summary$locus_seq %>% as.data.frame() %>% group_by_all() %>% count()
    summary(cs_count)
        #       .                   n         
        #  Length:164         Min.   : 2.000  
        #  Class :character   1st Qu.: 3.000  
        #  Mode  :character   Median : 6.000  
        #                     Mean   : 5.524  
        #                     3rd Qu.: 7.000  
        #                     Max.   :10.000  

# 5. get the CS with more than one SNP
    multivar_cs <- filter(cs_size_summary, num_variants>=2)
    multivar_count <- multivar_cs$locus_seq %>% as.data.frame() %>% group_by_all() %>% count()

# 6. Get the CS with only one SNP from the one with pip > .95
    singlevar_cs_95 <- final_results_loci %>% filter(locus_cs_id != multivar_cs$locus_cs_id) 
    nrow(singlevar_cs_95)

    singlevar_cs_95$SNP %>% unique() %>% length()

# 7. get the CS with only one SNP from the one with pip >.95 and lfsr <.01
    singlevar_cs_95_01 <- filtered_lfsr_df %>% filter(locus_cs_id != multivar_cs$locus_cs_id) 
    nrow(singlevar_cs_95_01)
    singlevar_cs_95_01$SNP %>% unique() %>% length()

# 8. get the CS with only one SNP from the one with pip >.95 and lfsr <.01 with at least one non Na traits (lfsr)
    singlevar_cs_95_01_noNA <- final_significant_results %>% filter(locus_cs_id != multivar_cs$locus_cs_id) 
    nrow(singlevar_cs_95_01_noNA)
    singlevar_cs_95_01_noNA$SNP %>% unique() %>% length()

# 9. get the number of non-NA CS in each trait
    as.data.frame(847- colSums(is.na(final_significant_results[,6:26])))
```

------------------------------------
[1]: <https://doi.org/10.1371/journal.pgen.1010299>
[2]: <https://github.com/stephenslab/mvsusieR/tree/master?tab=readme-ov-file>
[3]: <https://stephenslab.github.io/mvsusieR/articles/mvsusie_intro.html>
[4]: <wang.gao@columbia.edu>
