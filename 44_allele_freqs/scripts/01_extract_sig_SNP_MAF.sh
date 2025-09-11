#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --time=0:15:00
#SBATCH --mem=10G
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice
#SBATCH --partition=bigmem

cd /staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/44_allele_freqs/scripts

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

module purge
module load cluster/wice/bigmem
module load R/4.4.2-gfbf-2024a

eval "$(conda shell.bash hook)"
conda activate plink

INPUT_PLINK="../../16_pca_admixture_imputed/pca_with_htlv_only/htlv_no_fams"
SNP_DIR="../input/genome_suggestive_snps/"
OUTPUT_DIR="../output/snp_mafs/"
MAF_DIR="../output/low_maf/"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$MAF_DIR"

generate_maf_flag_table() {
    local frq_file=$1
    local raw_file=$2
    local out_file=$3

    Rscript Rscripts/generate_maf_table.R "$frq_file" "$raw_file" "$out_file"
}

# Main loop: process all SNP files
for SNP_FILE in "$SNP_DIR"/*.txt; do
    BASENAME=$(basename "$SNP_FILE" .txt)

    FORMATTED_SNP_FILE="$OUTPUT_DIR/${BASENAME}_formatted.txt"
    awk '{print $1 ":" $2}' "$SNP_FILE" > "$FORMATTED_SNP_FILE"

    OUTPUT_PLINK="$OUTPUT_DIR/$BASENAME"
    
    # Generate MAF (.frq) and genotype data (.raw)
    plink --bfile "$INPUT_PLINK" \
          --extract "$FORMATTED_SNP_FILE" \
          --freq \
          --recode A \
          --allow-no-sex \
          --out "$OUTPUT_PLINK"

    # File paths
    FRQ_FILE="$OUTPUT_PLINK.frq"
    RAW_FILE="$OUTPUT_PLINK.raw"
    MAF_OUTPUT_FILE="$MAF_DIR/${BASENAME}_maf_flags.tsv"

    # Call R script to produce final table
    generate_maf_flag_table "$FRQ_FILE" "$RAW_FILE" "$MAF_OUTPUT_FILE"

    echo "Processed $BASENAME → MAF + genotype table saved to $MAF_OUTPUT_FILE"
done

