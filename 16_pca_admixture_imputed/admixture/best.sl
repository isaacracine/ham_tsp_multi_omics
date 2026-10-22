#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCG –cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice

cd 16_pca_admixture_imputed/admixture

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
module purge
eval "$(conda shell.bash hook)"
conda activate admixture

bash admix_best_script --use-conda
