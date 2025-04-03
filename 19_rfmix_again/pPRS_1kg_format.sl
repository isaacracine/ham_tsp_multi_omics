#!/usr/bin/env bash

#SBATCH --nodes=2
#SBATCH --time=16:00:00
#SBATCH --mem=110G
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice
#SBATCH --ntasks-per-node=11


cd /staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/19_rfmix_again

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
module purge
eval "$(conda shell.bash hook)"
conda activate py-popgen


for i in {1..22}; 
do
bcftools view --threads 22 -O z -R ../24_pPRS/adding_LA/input/renamed_chr${i}_no_fam.vcf.gz -o ../24_pPRS/LA_1kg/input/1kg_chr${i}.vcf.gz 1kg_hg38/chr${i}_1kg38_filtered.recode.vcf.gz;
done
