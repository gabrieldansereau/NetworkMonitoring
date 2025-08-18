#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=01:15:00
#SBATCH --job-name=04_focal_array
#SBATCH --output=jobs/out/%x-%J.out

module load julia/1.11.3

cd $HOME/projects/def-tpoisot/gabdans/NetworkMonitoring

export PROGRESS_BARS_DT=60

julia --project -e 'const NREP = 50; const OUTDIR = "efficiency"; include("scripts/04_focal_sims.jl")'