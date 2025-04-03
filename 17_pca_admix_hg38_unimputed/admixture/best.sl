#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCG –cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=10000
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice

cd /staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/17_pca_admix_hg38_unimputed/admixture

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
module purge
eval "$(conda shell.bash hook)"
conda activate admixture

bash admix_best_script --use-conda
