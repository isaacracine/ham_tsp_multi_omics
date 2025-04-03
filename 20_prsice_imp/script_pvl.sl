#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mem=40000
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice

cd /staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/20_prsice_imp

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
module purge
eval "$(conda shell.bash hook)"
conda activate ibp_2022

Rscript PRSice.R \
    --prsice PRSice_linux \
    --base ss_no_hla.txt \
    --target no_amr_merged_all \
    --binary-target F \
    --pheno pheno_pvl.txt \
    --A1 ref \
    --A2 var \
    --seed 17 \
    --stat OR \
    --or \
    --bar-levels 5e-8,1e-5,0.001,0.005,0.01,0.05,0.1,0.5 \
    --fastscore T\
    --quantile 20 \
    --chr chr \
    --snp SNP \
    --bp loc_hg38 \
    --ld no_amr_merged_all \
    --clump-kb 250 \
    --clump-p 1.000000 \
    --clump-r2 0.100000 \
    --all-score \
    --print-snp \
    --no-regress \
    --out results_apr_11_no_amr_pvl_no_hla

#--print-snp \
# --cov hg38_unrelated/merged_covs.txt \
#--cov-col @PC_[1-4] \

