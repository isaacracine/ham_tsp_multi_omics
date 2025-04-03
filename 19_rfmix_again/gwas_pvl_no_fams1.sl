#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=70G
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice



cd /staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/19_rfmix_again

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
module purge
eval "$(conda shell.bash hook)"
conda activate py-popgen

for CHR in {1..6};
do
python Tractor/RunTractor.py \
      --hapdose tracts_4ancs_chr${CHR}_no_fam \
      --phe upd_covs_tractor/covs_with_age_tractor_pvl_no_fams.txt \
      --method linear \
      --out ss_4_anc_pvl_chr${CHR}_no_fam.tsv ;
done

