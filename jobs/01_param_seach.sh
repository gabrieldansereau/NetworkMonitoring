#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=20G
#SBATCH --time=01:00:00
#SBATCH --job-name=01
#SBATCH --output=jobs/out/%x-%J.out

module load julia/1.11.3

cd $HOME/projects/def-tpoisot/gabdans/NetworkMonitoring

srun hostname -s > hostfile
sleep 5
julia --project --machine-file ./hostfile ./scripts/01_param_search.jl
