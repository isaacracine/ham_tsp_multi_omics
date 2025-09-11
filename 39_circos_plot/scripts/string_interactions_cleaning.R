library(data.table)
library(dplyr)
library(circlize)
library(topr)

######### String interactions

##read in data

#read in the data
string_path <- "../input/string_interactions_short.tsv"
string_inters <- data.frame(fread(string_path, sep = "\t", header = TRUE))

#rename the cols
new_names <- c('node1')
colnames(string_inters)[1] <- new_names

#rename the protein ids
# Remove the prefix '9606.' using sub()
string_inters$node1_string_id <- sub("^9606\\.", "", string_inters$node1_string_id)
string_inters$node2_string_id <- sub("^9606\\.", "", string_inters$node2_string_id)

##write unique protein ids

# #extract unique genes to get chromosome positions from ensmeble online viewer
# # Combine the first two columns
# combined <- c(string_inters[[3]], string_inters[[4]])
# 
# # Get unique values
# unique_values <- unique(cleaned_values)
# 
# # Save to a text file
# id_out_path = "../output/unique_string_ensemble_ids.txt"
# writeLines(unique_values, id_out_path)


## read in ensembl results for chromosome locations from protein IDs
# read in the file
chr_pos_path = '../input/string_chrom_pos.csv'
chr_pos_full = read.csv(chr_pos_path, header = TRUE)

#keep chromosome, start, end and protein ID
chr_pos <- chr_pos_full[,c(8,3,4,5)]
colnames(chr_pos)[1]

#merge by protein ID: string_inters chr_pos
#merge first node info
si_ch <- merge(string_inters, chr_pos, by.x ='node1_string_id', by.y ='Protein.stable.ID', all.x = TRUE)
colnames(si_ch)[c(14:16)] <- c('chr_node1','gene_start_node1','gene_end_node1')
#merge second node info
si_ch <- merge(si_ch, chr_pos, by.x ='node2_string_id', by.y ='Protein.stable.ID', all.x = TRUE)
colnames(si_ch)[c(17:19)] <- c('chr_node2','gene_start_node2','gene_end_node2')

#find which proteins have NAs
node1_nas <- si_ch$node1_string_id[is.na(si_ch$chr_node1)]
node2_nas <- si_ch$node2_string_id[is.na(si_ch$chr_node2)]
miss_prot_id <- unique(c(node1_nas, node2_nas))

#manually add their information
#reference: uniprot -> then click on genomic location tab
mpd_1 <- c(miss_prot_id[1], 14, 54844020, 54902663) #14:54,844,020 - 54,902,663
mpd_2 <- c(miss_prot_id[2], 1, 32274330, 32285713) #1:32,274,330 - 32,285,713
mpd_3 <- c(miss_prot_id[3], 2, 268998, 271898) #2:268,998 - 271,898
mpd_4 <- c(miss_prot_id[4], 2, 233206, 262696) #2:233,206 - 262,696
missing_prot_info <- data.frame(rbind(mpd_1, mpd_2, mpd_3, mpd_4))
colnames(missing_prot_info) <- c('prot_id', 'chr','gene_start','gene_end')

#replace the missing information
for (i in 1:nrow(missing_prot_info)) {
  #if the prot name matches, replace chr, start and end
  si_ch$chr_node1[si_ch$node1_string_id==missing_prot_info$prot_id[i]] <- missing_prot_info$chr[i]
  si_ch$gene_start_node1[si_ch$node1_string_id==missing_prot_info$prot_id[i]] <- as.integer(missing_prot_info$gene_start[i])
  si_ch$gene_end_node1[si_ch$node1_string_id==missing_prot_info$prot_id[i]] <- as.integer(missing_prot_info$gene_end[i])
  
  si_ch$chr_node2[si_ch$node2_string_id==missing_prot_info$prot_id[i]] <- missing_prot_info$chr[i]
  si_ch$gene_start_node2[si_ch$node2_string_id==missing_prot_info$prot_id[i]] <- as.integer(missing_prot_info$gene_start[i])
  si_ch$gene_end_node2[si_ch$node2_string_id==missing_prot_info$prot_id[i]] <- as.integer(missing_prot_info$gene_end[i])
}

#make chr columns a string
si_ch$chr_node1 <- paste0("chr", si_ch$chr_node1)
si_ch$chr_node2 <- paste0("chr", si_ch$chr_node2)

##export cleaned data
cleaned_string_inter = '../input/string_interactions_cleaned.csv'
write.csv(si_ch, file = cleaned_string_inter, row.names = FALSE)


#### Now try plotting
# Initialize the circular plot
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))

# Add track for chromosomes
circos.trackPlotRegion(
  factors = unique(c(si_ch$chr_node1, si_ch$chr_node2)), 
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    #circos.text(CELL_META$xcenter, 1, chr, facing = "clockwise", adj = c(0, 0.5))
  }
)

# Add links with width proportional to the score
for (i in seq_len(nrow(si_ch))) {
  circos.link(
    sector.index1 = si_ch$chr_node1[i], point1 = c(si_ch$gene_start_node1[i], si_ch$gene_end_node1[i]),
    sector.index2 = si_ch$chr_node2[i], point2 = c(si_ch$gene_start_node2[i], si_ch$gene_end_node2[i]),
    col = rgb(1, 0, 0, alpha = 0.2), #si_ch$combined_score[i]), # Color transparency proportional to score
    lwd = si_ch$combined_score[i]               # Line width scaled by score
  )
}

# Clear the circular plot
circos.clear()


