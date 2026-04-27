#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --array=1-50
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=07_efficiency
#SBATCH --output=jobs/out/%x-%J.out

module load julia/1.11.3

cd $HOME/projects/def-tpoisot/gabdans/NetworkMonitoring

export PROGRESS_BARS_DT=60

julia --project -e 'const idn = 0; const NREP = 50; const OUTDIR = "efficiency"; include("scripts/04_focal_sims.jl")'
julia --project -e 'const idn = 50; const NREP = 50; const OUTDIR = "efficiency"; include("scripts/04_focal_sims.jl")'
julia --project -e 'const idn = 100; const NREP = 50; const OUTDIR = "efficiency"; include("scripts/04_focal_sims.jl")'
julia --project -e 'const idn = 150; const NREP = 50; const OUTDIR = "efficiency"; include("scripts/04_focal_sims.jl")'