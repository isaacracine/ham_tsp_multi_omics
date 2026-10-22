#this is the script for running a LassoSum analysis on the target
#database
#the install command needs to be run every new interactive session
library(devtools)
#install_github("tshmak/lassosum")

library(dplyr)
library(lassosum)
library(data.table)
library(methods)
library(magrittr)
library(parallel)

set.seed(17)

output_prefix = "apr_18_ham_ss_htlv"

#read in the files
sum.stat <- "input/sum_stats_hg38_pos.txt"
bfile <- "input/htlv_merged_with_fams_no_amr"

# Read in and process the covariates

#define the LD panel to use
ld.file <- "EUR.hg38"

# Read in the target phenotype file
target.pheno <- fread("input/pheno.txt")


#our names are messed up, rename
colnames(target.pheno) = c("FID", "IID", "pheno")


# Read in the summary statistics
ss <- fread(sum.stat)


#753 HAM/TSP patients and 899 asymptomatic 
ss$N <- 1652

#########
# Remove P-value = 0, which causes problem in the transformation
ss <- ss[!P == 0]



###################################################


#i Transform the P-values into correlation
cor <- p2cor(p = ss$P, #p_vals,#ORIGINAL: pull(ss, snakemake@params[[9]]),
        n = ss$N, 
        sign = ss$beta#log(as.numeric(ss$OR))
        )

#obtain fam file
fam <- fread(paste0(bfile, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]


bfile_ref <- "input/no_amr_merged_with_htlv_fams"
#ref_fam <- fread(paste0(bfile_ref, ".fam"), header = F)

#ref_fam[,6] <- NA
#write.table(ref_fam, paste0(bfile_ref, ".fam"), row.names = F, col.names = F, quote = F)

# Run the lassosum pipeline
# The cluster parameter is used for multi-threading
# You can ignore that if you do not wish to perform multi-threaded processing
out <- lassosum.pipeline( 
    cor = cor,
    chr = ss$chr,
    pos = ss$loc_hg38,
    A1 = ss$ref,
    A2 = ss$var,    
    ref.bfile = bfile_ref,
    test.bfile = bfile,
    LDblocks = ld.file,
    lambda=exp(seq(log(0.001), log(0.2), length.out=40)),
    s=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) 
    #cluster=cl
)
# we save the validation output
# we use validation rather than pseudovalidation, 
# because we want to allow for covariates 
target.res <- validate(out, pheno = as.data.frame(target.pheno), plot = FALSE)

#covar=as.data.frame(cov)

#save the best best s, lambda and pgs 
best_s <- target.res$best.s
best_lambda <- target.res$best.lambda
results_pgs <- target.res$results.table

#save the validation table from the results
validation_table <- target.res$validation.table

#save the validation plot
pdf(paste("output/", output_prefix,"_validation_plot.pdf", sep = '')) 
#Model validation refers to the process of confirming that the model actually achieves its intended purpose
plot(target.res)#, ylab = "Validation (correlation)")
# Close the pdf file
dev.off() 
#the best s and lambda will correspond to the highest validation
#value 


results_snp <- out$sumstats
results_beta <- target.res$best.beta
write.table(results_beta, quote = FALSE, row.names = FALSE, col.names = FALSE, paste("output/", output_prefix,"_betas.txt", sep = ''))
write.table(results_snp, quote = FALSE, row.names = FALSE, col.names = FALSE, paste("output/", output_prefix,"_snps.txt", sep = ''))

#save the result from the validation
write.table(best_s, quote = FALSE, row.names = FALSE, col.names = FALSE, paste("output/", output_prefix,"_best_s.txt", sep = ''))
write.table(best_lambda, quote = FALSE, row.names = FALSE, col.names = FALSE, paste("output/", output_prefix,"_best_lambda.txt", sep = ''))
write.table(results_pgs, quote = FALSE, row.names = FALSE,paste("output/", output_prefix,"_scores.txt", sep = ''))
write.table(validation_table, quote = FALSE, row.names = FALSE, paste("output/", output_prefix,"_validation_results.txt", sep = ''))


###obtain the prs for the best 4 combinations of the 
#thershold parameter, S, and shrinkage parameter, lambda
validation_table$Rank <- rank(-validation_table$value)
best_params <- validation_table[order(validation_table$Rank),]
write.table(best_params, quote = FALSE, row.names = FALSE, paste("output/", output_prefix,"_best_params.txt", sep = ''))


l1 <- best_params$lambda[1]
s1 <- best_params$s[1]
l2 <- best_params$lambda[2]
s2 <- best_params$s[2]
l3 <- best_params$lambda[3]
s3 <- best_params$s[3]
l4 <- best_params$lambda[4]
s4 <- best_params$s[4]

# Get the maximum R2
r2 <- max(target.res$validation.table$value)^2
write.table(r2, quote = FALSE, row.names = FALSE, col.names = FALSE, paste("output/", output_prefix,"_r_sqr_estimate.txt", sep = ''))
