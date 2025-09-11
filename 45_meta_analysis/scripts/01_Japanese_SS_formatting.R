library(dplyr)

#read in file
df <- read.table("/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/liftOver/gwas/sum_stats_hg38_pos.txt",
                 header = T)

#to calculate maf we need to consider the allele frequencies across
#cases and controls, as that is how the information is provided
#step 1: get sample size
df$N_case  <- with(df, CASE_HOMREF_CT + CASE_HET_CT + CASE_HOMVAR_CT)
df$N_ctrl  <- with(df, CTRL_HOMREF_CT + CTRL_HET_CT + CTRL_HOMVAR_CT)
df$N_all   <- df$N_case + df$N_ctrl

#step 2: case/control variant-allele frequencies
df$var_case_freq_calc <- with(df, (2*CASE_HOMVAR_CT + CASE_HET_CT) / (2*N_case))
df$var_ctrl_freq_calc <- with(df, (2*CTRL_HOMVAR_CT + CTRL_HET_CT) / (2*N_ctrl))

#step 3: get the overall MAF (since 'var' is the minor allele)
df$MAF_all <- with(df, {
  homvar_all <- CASE_HOMVAR_CT + CTRL_HOMVAR_CT
  het_all    <- CASE_HET_CT   + CTRL_HET_CT
  (2*homvar_all + het_all) / (2*N_all)
})



#### select the necessary columns 
#SNP     A1     A2     BETA     SE     P     N    FREQ
col_names <- c('SNP', 'ref', 'var', 
               'beta', 'ADD_SE', 'P', 'N_all', 'MAF_all')

# select the columns
ss <- df |> select(all_of(col_names))

# rename columns to match MR-MEGA's column names
#A1: effect allele 
new_names <- c("SNP", "A2", "A1", 
               "BETA", "SE",
               "P", "N", "FREQ")

colnames(ss) <- new_names


#save formatted SS
write.table(ss,
            "/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/45_meta_analysis/input/japanese_ss_hg38.txt",
            quote = F, sep = ' ', col.names = T, row.names = F)
