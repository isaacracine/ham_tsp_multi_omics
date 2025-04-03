#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --time=2:30:00
#SBATCH --mem=120G
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice

cd /staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/39_circos_plot/scripts

module load cluster/genius/batch
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

#Rscript df_cleaning.R
Rscript circos_plot.R 
