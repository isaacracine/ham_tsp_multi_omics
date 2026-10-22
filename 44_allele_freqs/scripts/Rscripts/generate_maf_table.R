args <- commandArgs(trailingOnly = TRUE)
frq_file <- args[1]
raw_file <- args[2]
out_file <- args[3]

# Read .frq file
frq <- read.table(frq_file, header = TRUE, stringsAsFactors = FALSE)
frq$SNP_ID <- frq$SNP

# Read .raw file (genotypes), prevent column name mangling
raw <- read.table(raw_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
sample_ids <- raw$IID
geno <- raw[, -(1:6)]  # remove FID, IID, PAT, MAT, SEX, PHENOTYPE

# Strip allele suffix from column names
raw_snp_names <- colnames(geno)
stripped_snp_names <- sub("_[ACGT]+$", "", raw_snp_names)

# Keep SNPs present in both files
keep_idx <- stripped_snp_names %in% frq$SNP
geno <- geno[, keep_idx, drop = FALSE]
stripped_snp_names <- stripped_snp_names[keep_idx]
frq_subset <- frq[match(stripped_snp_names, frq$SNP), ]

# Rename genotype columns to match index
colnames(geno) <- paste0("SNP_", seq_len(ncol(geno)))

# Initialize sample-wise output table
sample_matrix <- data.frame(sample_ID = sample_ids, stringsAsFactors = FALSE)

for (i in seq_along(colnames(geno))) {
  snp_name <- frq_subset$SNP_ID[i]
  maf <- frq_subset$MAF[i]
  dosage <- geno[[i]]

  a1 <- dosage
  a2 <- 2 - dosage

  # Create generic column block: SNP_X, MAF_X, A1_X, A2_X
  block <- data.frame(
    snp_col = rep(snp_name, length(dosage)),
    maf_col = rep(maf, length(dosage)),
    A1 = a1,
    A2 = a2
  )

  names(block) <- c(
    paste0("SNP_", i),
    paste0("MAF_", i),
    paste0("A1_", i),
    paste0("A2_", i)
  )

  sample_matrix <- cbind(sample_matrix, block)
}

# Write to output
write.table(sample_matrix, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

