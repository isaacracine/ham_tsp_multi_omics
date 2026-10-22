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
cd 45_meta_analysis/

# Set metal path as variables
METAL=45_meta_analysis/generic-metal/metal
METAL_SCRIPT=45_meta_analysis/scripts/metal_script.txt

# Run metal script
$METAL $METAL_SCRIPT
