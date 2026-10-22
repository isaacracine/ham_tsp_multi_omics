
ancs <- c("afr", "amr", "asn", "eur")
phenos <- c("ham", "pvl")

for(i in ancs){
	for(j in phenos){
		data <- read.table(paste('../', i, "_", j, "_", 'ss.txt', sep = ''), header = TRUE)
		data <- data[!is.na(data$P),]
		data$chisq <- qchisq(1 - data$P, df =1)
		lambda_gc <- median(data$chisq) / qchisq(0.5, df = 1)
		cat("Genomic inflation factor (λGC):", lambda_gc, "\n", file = paste("../", i, "_", j, "_genomic_inflation_factor.txt", sep =""))

	}
}


