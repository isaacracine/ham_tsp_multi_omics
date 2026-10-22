#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --mem=100G
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice



cd 19_rfmix

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
module purge
eval "$(conda shell.bash hook)"
conda activate py-popgen

bash rfmix_script --use-conda
