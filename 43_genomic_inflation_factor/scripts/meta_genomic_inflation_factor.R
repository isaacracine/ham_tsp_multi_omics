data <- read.table('/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/45_meta_analysis/output/meta_ham_ss.txt', header = TRUE)
data <- data[!is.na(data$P),]
data$chisq <- qchisq(1 - data$P, df =1)
lambda_gc <- median(data$chisq) / qchisq(0.5, df = 1)

cat("Genomic inflation factor (λGC):", lambda_gc, "\n", file = "/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/43_genomic_inflation_factor/output/meta_ham_genomic_inflation_factor.txt")

