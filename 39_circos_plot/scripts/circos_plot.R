library(data.table)
library(dplyr)
library(tidyr)
library(circlize)
library(stringr)
library(ComplexHeatmap)

######### GWAS

#read in the data
regenie_file_path <- "../../35_regenie/output/ham_ss_regenie.txt"
r_ss_df <- data.frame(fread(regenie_file_path, sep = " ", header = TRUE))
r_drops <- c("REF", "ALT", "SNP","N", "BETA","SE")
#downsample for testing
r_ss <- r_ss_df[ , !(names(r_ss_df) %in% r_drops)]
#r_ss <- r_ss[seq(1, nrow(r_ss), by = 10000), ]

tractor_file_path <- "../../19_rfmix_again/ham_ss_tractor.txt"
t_ss_df <- data.frame(fread(tractor_file_path, sep = " ", header = TRUE))
t_drops <- c("SNP", "REF", "ALT", "BETA_AFR", "BETA_AMR", "BETA_ASN", "BETA_EUR")
#downsample for testing
t_ss <- t_ss_df[ , !(names(t_ss_df) %in% t_drops)]
#t_ss <- t_ss[seq(1, nrow(t_ss), by = 10000), ]

#format into one data frame
#join two data frames on columns CHR and POS
tot_ss <- merge(r_ss, t_ss, by = c("CHR", "POS"), all = TRUE)

print('done GWAS')

############ STRING interactions
si_path = '../input/string_interactions_cleaned.csv'
string_int = read.csv(si_path)

# Combine both columns into one and count occurrences
gene_counts <- string_int %>%
  select(node1_string_id, node2_string_id) %>%
  pivot_longer(cols = everything(), values_to = "ID") %>%
  count(ID, name = "Count")

# Calculate the top 10% most connected genes
hub_genes <- gene_counts %>%
  arrange(desc(Count)) %>%
  slice(1:ceiling(n() * 0.2))

# Filter rows where the gene is in col1 or col2 and extract locations
locations <- string_int %>%
  filter(node1_string_id %in% hub_genes$ID | node2_string_id %in% hub_genes$ID) %>%
  mutate(
    gene = case_when(
      node1_string_id %in% hub_genes$ID ~ node1_string_id,
      node2_string_id %in% hub_genes$ID ~ node2_string_id
    ),
    node = case_when(
      gene == node1_string_id ~ node1,
      gene == node2_string_id ~ node2
    ),
    chr = case_when(
      gene == node1_string_id ~ chr_node1,
      gene == node2_string_id ~ chr_node2
    ),
    gene_start = case_when(
      gene == node1_string_id ~ gene_start_node1,
      gene == node2_string_id ~ gene_start_node2
    ),
    gene_end = case_when(
      gene == node1_string_id ~ gene_end_node1,
      gene == node2_string_id ~ gene_end_node2
    )
  ) %>%
  select(gene, node, chr, gene_start, gene_end) %>%
  unique()

# get the information for AIM2
aim2_loc <-unique(string_int[string_int$node1=='AIM2',c('node1_string_id', 'node1', 'chr_node1', 'gene_start_node1', 'gene_end_node1')])
colnames(aim2_loc) <- colnames(locations)

#combine with hub locations 
locations <- rbind(locations, aim2_loc)

#reorganize to be in bed format
locations <- locations[c(3,4,5,2,1)]

print('done STRING')

############ Trasncriptomics results
tr_res_path = '../input/transcript_results_cleaned.csv'
tr_res = read.csv(tr_res_path)
tr_res$CHR <- paste("chr", tr_res$CHR, sep="")

#remove x and y chrom entries
tr_res <- tr_res[!(tr_res$CHR %in% c("chrX", "chrY")), ]


print('done TRANSCRIPTOMICS')

################## Circos Official Plot
##### Further cleaning
gwas_data <- tot_ss

# Add 'chr' prefix to the chromosome names
gwas_data$CHR <- paste0("chr", gwas_data$CHR)

# Calculate chromosome ranges
#gwas
chr_ranges_gwas <- aggregate(POS ~ CHR, gwas_data, function(x) c(min(x), max(x)))
chr_ranges_gwas <- do.call(data.frame, chr_ranges_gwas)
names(chr_ranges_gwas) <- c("chromosome", "start", "end")

#transcriptomics
chr_ranges_trans <- aggregate(
  tr_res[, c("CHRLOC", "CHRLOCEND")], 
  by = list(CHR = tr_res$CHR), 
  FUN = function(x) if (identical(x, tr_res$CHRLOC)) min(x) else max(x)
)

#set bar width to use in transcriptomics and add it to CHRLOCEND
bar_width <- 1e7
chr_ranges_trans$CHRLOCEND <- chr_ranges_trans$CHRLOCEND + bar_width + 100 #add 100 to make sure enough space
chr_ranges_trans <- do.call(data.frame, chr_ranges_trans)
names(chr_ranges_trans) <- c("chromosome", "start", "end")

#select min/max from both gwas and transcriptomics 
combined_df <- rbind(
  chr_ranges_gwas[, c("chromosome", "start", "end")], 
  chr_ranges_trans[, c("chromosome", "start", "end")]
)

result_df <- combined_df %>%
  group_by(chromosome) %>%
  summarise(
    start = min(start),
    end = max(end)
  ) %>%
  ungroup()

# Convert to a data frame if needed
chr_ranges <- as.data.frame(result_df)
colnames(chr_ranges) <- c("chromosome","start","end")

#extract the chr number of each and sort numerially
chr_num <- list()
for (i in 1:nrow(chr_ranges)) {
  chr_num[[length(chr_num)+1]] <- as.integer(str_split(chr_ranges$chromosome, "chr")[[i]][2])
  
}
chr_ranges$chr_num <- unlist(chr_num)
chr_ranges <- chr_ranges[order(chr_ranges$chr_num),] 


#merge with tr_res and remove out of box values
# Merge df2 with chr_ranges based on Chromosome
merged_df <- merge(tr_res, chr_ranges, by.x = 'CHR', by.y = "chromosome")

# Filter out rows where:
# - Start is smaller than the corresponding chromosome start
# - End is larger than the corresponding chromosome end
# - Significant for p-value 
filtered_df <- merged_df %>%
  filter(CHRLOC >= start & CHRLOCEND <= end) %>%
  filter(P.Value < 0.05)

#extend front of chromosome tracks
chr_ranges$start <- chr_ranges$start - bar_width - 100 #add 100 to make sure enough space


#make it so the chromsome tracks plot numerically
chr_ranges$chromosome <- factor(
  chr_ranges$chromosome, 
  levels = paste('chr',chr_ranges$chr_num, sep = '')
)

########Begin acutal plotting

# Adjust Circos Parameters for compact display
circos.par(
  gap.degree = 1.5,  # Reduce gaps between chromosomes to make more space for tracks
  track.margin = c(0.01, 0.01),  # Adjust margins between tracks (smaller values for less space)
  cell.padding = c(0, 0, 0, 0),  # Reduce cell padding to fit more in the plot
  start.degree = 90  # Adjust starting degree for better orientation
)

#Initialize png
#Open PNG device, setting the dimensions to 800x800
png("../output/circos_plot.png", width = 1600, height = 1600)


# Initialize Circos Plot
#circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
# circos.initialize(factors = chr_ranges$chromosome,
#                   xlim = cbind(chr_ranges$start, chr_ranges$end))
circos.genomicInitialize(chr_ranges, plotType = NULL)

#plot chromosome names
circos.track(
  track.height = 0.0001,
  ylim = c(0, 1),  # Define the y-limits for the track
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,  # Get the center of each chromosome sector
      CELL_META$ylim[2] + mm_y(-15),  # Position the text above the track
      gsub(".*chr", "", CELL_META$sector.index),  # Extract the chromosome name
      cex = 2,  # Font size
      niceFacing = TRUE  # Ensure the text is readable and faces outward
    )
  },
  cell.padding = c(0, 0, 0, 0),  # No padding between the track and the plot
  bg.border = NA  # No border around the track
)

#plot hub genes
circos.genomicLabels(locations, 
                     labels.column = 4,  # Column containing the label text
                     side = "outside",   # Place labels outside the circle
                     col = "blue",       # Label color
                     line_col = "gray",  # Connector line color
                     cex = 1.5)          # Font size for the labels

circos.par(
  track.margin = c(-0.01, -0.01)  # Adjust margins between tracks (smaller values for less space)
)

# Extend the connector lines across tracks
for (i in 1:nrow(locations)) {
  circos.link(
    sector.index1 = locations$chr[i],
    point1 = c(locations$gene_start[i], locations$gene_end[i]),
    sector.index2 = locations$chr[i],
    point2 = c(locations$gene_start[i], locations$gene_end[i]),
    col = "gray",      # Link color
    lwd = 0.7,           # Line width
    border = "gray",   # Border for the link
    h = 0.46 # set the height of the line
  )
}

circos.par(
  track.margin = c(0.01, 0.01)  # Adjust margins between tracks (smaller values for less space)
)

#set color for transcript pval
#set the color map
col_fun <- colorRamp2(c(0.05, 1e-5, 5e-8), c("blue", "yellow", "red"), transparency = 0.25)

# transcriptomics bar chart 
circos.track(
  ylim = range(min(filtered_df$logFC)-0.35, max(filtered_df$logFC)+0.45), 
  track.height = 0.2, # Adjust track height as needed
  panel.fun = function(region, value, ...) {
    chr = get.current.chromosome()
    # Subset data for the current chromosome
    chr_data = filtered_df[filtered_df$CHR == chr, ]
    # Use circos.barplot to add bars
    circos.barplot(
      -chr_data$logFC, 
      chr_data$CHRLOC, 
      bar_width = bar_width,
      col = col_fun(chr_data$P.Value), # Color based on P.val
      border = NA
    )
    # Add a horizontal line at 0
    circos.lines(c(min(chr_data$CHRLOC), max(chr_data$CHRLOC)), 
                 c(0, 0), 
                 col = "black", 
                 lty = 1, # Solid line
                 lwd = 0.5) # Line width
  }
)

# Now create a separate continuous color legend using ComplexHeatmap's Legend
lgd = Legend(
  col_fun = col_fun, 
  title = "p-value color scale", 
  at = c(5e-8, 0.05),
  #new
  grid_height = unit(15, "mm"), # Adjust the size of legend items
  grid_width = unit(10, "mm"),
  title_gp = gpar(fontsize = 20), # Adjust the title size
  labels_gp = gpar(fontsize = 18)) # Adjust the labels size

# Draw the legend
draw(lgd, x = unit(0.9, "npc"), y = unit(0.1, "npc"))  # Use grid units


# Function to create Manhattan plot track
create_manhattan_track <- function(gwas_data, p_value_column, track_color, threshold = -log10(5e-8), sug_threshold = -log10(1e-5), track_height = 0.2) {
  
  # Ensure the column exists in the data
  if (!(p_value_column %in% colnames(gwas_data))) {
    stop(paste("Column", p_value_column, "does not exist in the data"))
  }
  
  # Transform p-values to -log10 scale, and remove NAs
  gwas_data$log_p <- -log10(gwas_data[[p_value_column]])
  
  # Remove NA values from the log_p column
  gwas_data <- gwas_data[!is.na(gwas_data$log_p), ]
  
  # Add p-value track
  circos.trackPlotRegion(
    factors = as.factor(gwas_data$CHR), 
    y = gwas_data$log_p,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      xlim <- CELL_META$xlim
      ylim <- CELL_META$ylim
    },
    ylim = c(0, max(gwas_data$log_p)),
    track.height = track_height
    #remove the border
    #bg.border = NA  # Remove default track border
  )
  
  # Add Points for SNPs with specified color
  circos.trackPoints(factors = as.factor(gwas_data$CHR), 
                     x = gwas_data$POS, 
                     y = gwas_data$log_p, 
                     pch = 16, col = track_color, cex = 0.5)
  
  # Add Significance Threshold Line
  # Loop through each chromosome in chr_ranges
  for (i in 1:nrow(chr_ranges)) {
    # Extract chromosome name, start, and end positions
    chr <- chr_ranges$chromosome[i]
    start <- chr_ranges$start[i]
    end <- chr_ranges$end[i]
    
    # Draw the significance threshold line
    circos.lines(
      x = c(start, end),
      y = c(threshold, threshold),
      col = "red",
      lwd = 0.5,
      sector.index = chr       # Draw in the correct sector
      #track.index = 3           # Explicitly specify track 3
    )
    
    # Draw the significance threshold line
    circos.lines(
      x = c(start, end),
      y = c(sug_threshold, sug_threshold),
      col = "blue",
      lwd = 0.5,
      sector.index = chr       # Draw in the correct sector
      #track.index = 3           # Explicitly specify track 3
    )
  }
    
}

#regnie
#gwas_data[1:500,'P'] <- 2e-10 # make it significant
create_manhattan_track(gwas_data, p_value_column = "P", track_color = "darkorange")

# #afr
# create_manhattan_track(gwas_data, p_value_column = "P_AFR", track_color = "red")
# 
# #amr
# create_manhattan_track(gwas_data, p_value_column = "P_AMR", track_color = "pink")
# 
# #asn
# create_manhattan_track(gwas_data, p_value_column = "P_ASN", track_color = "green")
# 
# #eur
# create_manhattan_track(gwas_data, p_value_column = "P_EUR", track_color = "black")



# Add links with width proportional to the score
for (i in seq_len(nrow(string_int))) {
  circos.link(
    sector.index1 = string_int$chr_node1[i], point1 = c(string_int$gene_start_node1[i], string_int$gene_end_node1[i]),
    sector.index2 = string_int$chr_node2[i], point2 = c(string_int$gene_start_node2[i], string_int$gene_end_node2[i]),
    #constanct color
    #col = rgb(1, 0, 0, alpha = 0.2), #si_ch$combined_score[i]), # Color transparency proportional to score
    #random color
    col = rand_color(nrow(string_int), transparency = 0.2),
    lwd = string_int$combined_score[i] *2               # Line width scaled by score
  )
}

# Clear Circos Plot when done
circos.clear()

#End plot
dev.off()

