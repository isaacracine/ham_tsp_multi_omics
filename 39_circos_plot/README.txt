Input:
* string_interactions_short.tsv: results from running STRING interactions with FUMA results from GWAS
* HAM_GWAS_transcriptomics_GSE29333.top.xlsx: results from transcriptomics analysis but multi-sheet excel file
* transcript_results.csv: the sheet entitled `Tattermusch2_GSE29333.top.table` was extracted
  * pip install csvkit 
  * in2csv your_file.xlsx --sheet "Tattermusch2_GSE29333.top.table" > transcript_results.csv
* ham_ss_regenie result were also used for creation of the plot
  
  
Scripts:
*circos_plot.R: plot used to construct the circos plot after formating the data
*transcriptomic_cleaning.R: script used to clean transcriptomic reads and for mapping the reads to their corresponding genome position
*string_interactions_cleaning.R: used to clean results from string analysis and to obatain their genomic positions


Processes:
* The output file labeled `unique_string_ensemble_ids.txt` were passed to the ensemble database to get the chromosome and positions. 
  * Ensemble Genes 113
  * Human genes (GRCh38.p14)
  * Within filters tab:  Input external references ID list [Max 500 advised] and select Protein stable ID(s) [e.g. ENSP00000000233] and pass in IDs from file
  * Attributes: add to get the chromsome and position, and protein ID so we can match them back up
  * 50 unique protein IDs were in the file, but only 46 IDed, will need to search for the remainders
  
  

