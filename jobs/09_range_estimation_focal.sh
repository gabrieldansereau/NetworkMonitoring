#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=02:30:00
#SBATCH --job-name=09_range_estimation_focal
#SBATCH --output=jobs/out/%x-%J.out

module load julia/1.11.3

cd $HOME/projects/def-tpoisot/gabdans/NetworkMonitoring

export PROGRESS_BARS_DT=60

julia --project -e 'const NREP = 100; const OUTDIR = "sim-range-estimation"; include("scripts/04_focal_sims.jl")'