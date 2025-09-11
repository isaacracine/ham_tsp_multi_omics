library(dplyr)

## Test for HAM

ham_file_list <- c('../output/pheno_low_maf/afr_ham_genom_sug_ss_maf_flags_pheno.txt',
                '../output/pheno_low_maf/eur_ham_genom_sug_ss_maf_flags_pheno.txt',
                '../output/pheno_low_maf/regenie_ham_genom_sug_ss_maf_flags_pheno.txt')

maf_thresh = 0.05

fischer_p_vals <- lapply(ham_file_list, function(file_names){
  # Read table
  t <- read.table(file_names,
                  header = T)
  
  # Get MAF cols
  t_maf <- t[, grepl("^MA", names(t))][1,] |> unlist()

  # Get the MAF cols with a low maf
  maf_names <- which(lapply(t_maf, function(x) x < maf_thresh) |>
                       unlist()) |> names() 
  
  # Find the position of column "i"
  col_index <- which(names(t) %in% maf_names)
  
  # From the column index, get the SNP name (1 to the left)
  # and the two alleles (1 and 2 cols to the right)
  # also include the phenotypes
  result_list <- lapply(col_index, function(col_name) {
    start <- max(col_name - 1, 1)
    end <- min(col_name + 2, ncol(t))
    t[, c(1:3,start:end)]
  })
  
  
  # Now we need to add a flag
  result_list2 <- lapply(result_list, function(df) {
    df$encoding <- as.factor(ifelse(df[[6]] == 0, "homo_ref",
                          ifelse(df[[6]] == 1, 'hetero', 'homo_alt')))
    return(df)
  })
  
  # Now we can run the fisher exact tests
  # Could also do this as a one sided test but must now if
  # protective or risk variant from sum stats (OR != 1)
  result_list3 <- lapply(result_list2, function(df){
    res <- table(df$encoding, df$ham) |> as.matrix() |> fisher.test()
    return(res$p.value)
  })
  
  # Now need to add the SNP name and the p-value to return
  snps <- lapply(result_list, function(res)
    {
      unname(unique(res[4]))[[1]]
  })

  # Calculate the FDR
  snp_fdr <- p.adjust(unlist(result_list3), method = "fdr")
  
  # Return as df
  df <- data.frame(
    snp = unlist(snps),
    p.val = unlist(result_list3),
    fdr = unlist(snp_fdr)
  )

  
  return(df)
})


fischer_p_vals

# Now need to save the list of dataframes as different files!
# Only the third has low enough MAFs
write.table(fischer_p_vals[[3]],
            "../output/low_maf_fischer/regnie_ham_genom_sug_ss_maf_fischer.csv",
            quote = F, sep = ',', col.names = T, row.names = F)


######## Continue testing for PVL

pvl_file_list <- c('../output/pheno_low_maf/afr_pvl_genom_sug_ss_maf_flags_pheno.txt',
                   '../output/pheno_low_maf/amr_pvl_genom_sug_ss_maf_flags_pheno.txt',
                   '../output/pheno_low_maf/eur_pvl_genom_sug_ss_maf_flags_pheno.txt',
                   '../output/pheno_low_maf/regenie_pvl_genom_sug_ss_maf_flags_pheno.txt')

anova_p_vals <- lapply(pvl_file_list, function(file_names){
  # Read table
  t <- read.table(file_names,
                  header = T)

  # drop rows with NA for pvl
  t <- t |> filter(!is.na(pvl))
  
  
  # log transform pvl if not already
  t$pvl <- log10(t$pvl)
  
  # Get MAF cols
  t_maf <- t[, grepl("^MA", names(t))][1,] |> unlist()
  
  # Get the MAF cols with a low maf
  maf_names <- which(lapply(t_maf, function(x) x < maf_thresh) |>
                       unlist()) |> names() 
  
  # Find the position of column "i"
  col_index <- which(names(t) %in% maf_names)
  
  # From the column index, get the SNP name (1 to the left)
  # and the two alleles (1 and 2 cols to the right)
  # also include the phenotypes
  result_list <- lapply(col_index, function(col_name) {
    start <- max(col_name - 1, 1)
    end <- min(col_name + 2, ncol(t))
    t[, c(1:3,start:end)]
  })
  
  
  # format the data to run ANOVA
  # needs to be formated so that genotype (1,2,3) is one column and the pvl is another!!
  # should again get the columns for each maf + allele (as above) and then just encode using
  # the dosage and the pvl
  
  # Now we need to add a flag
  result_list2 <- lapply(result_list, function(df) {
    df$encoding <- as.factor(ifelse(df[[6]] == 0, "homo_ref",
                                    ifelse(df[[6]] == 1, 'hetero', 'homo_alt')))
    return(df)
  })
  
  # Now can run ANOVA
  result_list3 <- lapply(result_list2, function(df)
  {
    res <- aov(df$pvl ~ factor(df$encoding))
    res <- cor.test(df$pvl, as.numeric(factor(df$encoding))-1, method = "spearman", exact = FALSE)
    return(res$p.value)
      #return(unlist(summary(res))[9])
  })

  # Now need to add the SNP name and the p-value to return
  snps <- lapply(result_list, function(res)
  {
    unname(unique(res[4]))[[1]]
  })
  
  # Calculate the FDR
  snp_fdr <- p.adjust(unlist(result_list3), method = "fdr")
  
  # Return as df
  df <- data.frame(
    snp = unlist(snps),
    p.val = unlist(result_list3),
    fdr = unlist(snp_fdr)
  )
  
})

anova_p_vals

#save results to a file
#afr_pvl
write.table(anova_p_vals[[1]],
            "../output/low_maf_anova/afr_pvl_genom_sug_ss_maf_flags_anova.csv",
            quote = F, sep = ',', col.names = T, row.names = F)

#amr_pvl -> nothing remained

#eur_pvl
write.table(anova_p_vals[[3]],
            "../output/low_maf_anova/eur_pvl_genom_sug_ss_maf_flags_anova.csv",
            quote = F, sep = ',', col.names = T, row.names = F)

#regenie_pvl
write.table(anova_p_vals[[4]],
            "../output/low_maf_anova/regenie_pvl_genom_sug_ss_maf_flags_anova.csv",
            quote = F, sep = ',', col.names = T, row.names = F)





