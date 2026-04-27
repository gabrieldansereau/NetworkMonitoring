#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --job-name=01_param_search
#SBATCH --output=jobs/out/%x-%J.out

module load julia/1.11.3

cd $HOME/projects/def-tpoisot/gabdans/NetworkMonitoring

export PROGRESS_BARS_DT=60

julia --project -e 'const NREP = 50; const OUTDIR = "sim-prop"; include("scripts/01_param_search.jl")'
