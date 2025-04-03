#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --mem=100G
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice



cd /staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/19_rfmix_again

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
module purge
eval "$(conda shell.bash hook)"
conda activate py-popgen

for CHR in {1..6};
do
python Tractor/ExtractTracts.py \
      --msp full_rf_ancsetry_${CHR}_no_fam \
      --vcf htlv/renamed_chr${CHR}_no_fam \
      --zipped \
      --output-path tracts_4ancs_chr${CHR}_no_fam \
      --num-ancs 4 ;
done
