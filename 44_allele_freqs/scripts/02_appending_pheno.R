
#load packages
library(tidyverse)

#set paths
base_path <- "/staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/44_allele_freqs/output/low_maf/"
ham_path <- "/staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/20_prsice_imp/pheno.txt"
pvl_path <- "/staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/20_prsice_imp/pheno_pvl.txt"
ss_list <-c("afr_ham_genom_sug_ss_maf_flags",
	"afr_pvl_genom_sug_ss_maf_flags",
	"amr_pvl_genom_sug_ss_maf_flags",
	"eur_ham_genom_sug_ss_maf_flags",
	"eur_pvl_genom_sug_ss_maf_flags",
	"regenie_ham_genom_sug_ss_maf_flags",
	"regenie_pvl_genom_sug_ss_maf_flags")

#read in phenos
ham <- read.table(ham_path, header = F, sep = ' ') 
pvl <- read.table(pvl_path, header = F, sep = ' ') 

#drop first column
ham <- ham[, -1]
pvl <- pvl[, -1]

#rename phenos to match
colnames(ham) <- c('sample_ID', 'ham')
colnames(pvl) <- c('sample_ID', 'pvl')

#do for loop to open sum stats, append and save
for(ss in ss_list){
	sum_stats <- read.table(paste(base_path, ss, ".tsv", sep = ''), header = T, sep = '\t') 
	sum_stats1 <- merge(pvl, sum_stats, by = 'sample_ID')
	sum_stats2 <- merge(ham, sum_stats1, by = 'sample_ID')
	write.table(sum_stats2, paste("/staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/44_allele_freqs/output/pheno_low_maf/", ss, "_pheno.txt", sep = ""), sep = " ", quote = F, row.names = F, col.names = T)
}




