library(dplyr)

#read in the file
df <- read.table("/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/35_regenie/output/ham_ss_regenie.txt",
                 header = T)

#convert beta to OR
df <- df |>   mutate(
  OR      = exp(BETA),
  OR_95L  = exp(BETA - 1.96 * SE),
  OR_95U  = exp(BETA + 1.96 * SE)
)

### need to get MAF!
# Read PLINK freq output
maf_df <- read.table("/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/45_meta_analysis/input/regenie_ham_maf.frq",
                     header = T)

# Make sure SNP IDs match between PLINK and your dataframe
# Suppose your GWAS summary stats are in `gwas_df` with column "MARKERNAME"
merged_df <- merge(df, maf_df |> select(SNP, MAF), by = "SNP")


### renaming
#SNP     A1     A2     BETA     SE     P     N    FREQ
#A1: effect allele 
col_names <- c('SNP', 'ALT', 'REF', 
               'BETA', 'SE',
               'P', 'N', 'MAF')

# select the columns
ss <- merged_df |> select(all_of(col_names))


# rename columns to match MR-MEGA's column names
new_names <- c("SNP", "A1", "A2", 
               "BETA", "SE",
               "P", "N", "FREQ")

colnames(ss) <- new_names

#save formatted SS
write.table(ss,
            "/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/45_meta_analysis/input/ham_regenie_ss_hg38.txt",
            quote = F, sep = ' ', col.names = T, row.names = F)

