#!/usr/bin/env bash

#SBATCH --job-name=mama_job
#SBATCH --nodes=1
#SBATCH --time=15:00:00
#SBATCH --mem=30G
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice
#SBATCH -o bash_logs/%x_%j.out
#SBATCH -e bash_logs/%x_%j.err

# Set path
cd /staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/45_meta_analysis/

# Set metal path as variables
METAL=/staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/45_meta_analysis/generic-metal/metal
METAL_SCRIPT=/staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/45_meta_analysis/scripts/metal_script.txt

# Run metal script
$METAL $METAL_SCRIPT
