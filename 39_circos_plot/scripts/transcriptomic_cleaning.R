library(data.table)
library(dplyr)
library(circlize)
library(illuminaHumanv4.db)

######### Transcriptomic data

#read in the data
transcript_path <- "../input/transcript_results.csv"
tr_res <- data.frame(fread(transcript_path, sep = ",", header = TRUE))

###see if can get gene location from illumina ID
# Function to annotate a dataframe column of Illumina IDs
annotate_probe_ids <- function(df, column_name) {
  # Extract probe IDs from the specified column
  probe_ids <- as.character(df[[column_name]])
  
  # Fetch annotations for the probe IDs
  annotations <- select(illuminaHumanv4.db, 
                        keys = probe_ids, 
                        keytype = "PROBEID", 
                        columns = c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND"))
  
  # Merge annotations back into the original dataframe
  merged_df <- merge(df, annotations, by.x = column_name, by.y = "PROBEID", all.x = TRUE)
  
  # Return the annotated dataframe
  return(merged_df)
}

# Annotate the dataframe
annotated_df <- annotate_probe_ids(tr_res, "ID")

# filter for illumina IDs where all entries have NA for chromosome
result <- annotated_df %>%
  group_by(ID) %>%
  filter(all(is.na(CHR))) %>%
  ungroup() %>%
  distinct(ID) %>%
  pull(ID)

#get the the IDs with no chromsome information
result <- annotated_df %>%
  filter(ID %in% result) 

#check if all have adj p-val = 1
nrow(result)
sum(result$P.Val<=1)

#YES all the p-values for these IDs without a chrom
#are equal to 1, so they can be ignored in th plot

###### Select the chromsome start and end to use 
#some of the position are negative -> antisense strad
#make them all absolute values
annotated_df$CHRLOC <- abs(annotated_df$CHRLOC)
annotated_df$CHRLOCEND <- abs(annotated_df$CHRLOCEND)


#do some filtering where we only take entries with:
  # not NA CHR
  # take the first enrtry of a unique ID
tr_res_chr_pos <- annotated_df[c(1,2,3,6,11,12,13,14)] %>%
  # Filter out rows where 'chr' or 'p-val' are NA
  filter(!is.na(CHR),!is.na(CHRLOC),!is.na(CHRLOCEND)) %>%
  # Remove duplicates while keeping the first occurrence of each 'ID'
  distinct(ID, .keep_all = TRUE)


######save the dataframe
id_out_path = '../input/transcript_results_cleaned.csv'
write.csv(tr_res_chr_pos, id_out_path)
