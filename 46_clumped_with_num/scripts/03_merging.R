# Read in files
args <- commandArgs(trailingOnly = TRUE)

clumped_ss <- read.table(args[1], header = TRUE, stringsAsFactors = FALSE)
counts <- read.table(args[2], header = FALSE, stringsAsFactors = FALSE)
output_file <- args[3]

# Add column names
colnames(counts) <- c("SNP", "CLUMP_COUNT")

# Merge
clump_merg <- merge(clumped_ss, counts, by = 'SNP')

# Reorder
# Make SNP the 3rd column
if(ncol(clump_merg) > 9)
{
  clump_merg <- clump_merg[, c(2, 3, 1, 4, 5, 6, 7, 8, 9, 10)]
}else
{
  clump_merg <- clump_merg[, c(2, 3, 1, 4, 5, 6, 7, 8)]
}

# Sort ascendingly by the p-value
clump_merg <- clump_merg[order(clump_merg$P), ]

# Write results
write.table(clump_merg, output_file, row.names = F, col.names = T, quote = F)

