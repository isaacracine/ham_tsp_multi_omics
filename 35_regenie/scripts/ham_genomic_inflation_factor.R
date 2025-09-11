

data <- read.table('../output/ham_ss_regenie.txt', header = TRUE)
data <- data[!is.na(data$P),]
data$chisq <- qchisq(1 - data$P, df =1)
lambda_gc <- median(data$chisq) / qchisq(0.5, df = 1)

cat("Genomic inflation factor (λGC):", lambda_gc, "\n", file = "../output/ham_regenie_genomic_inflation_factor.txt")



