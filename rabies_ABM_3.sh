#! /bin/bash
#SBATCH --array=1001-1125
#SBATCH --time=02:30:00
#SBATCH --cpus-per-task=1 
#SBATCH --account=ctb-tpoisot
#SBATCH --output=slurm\%x-%a.out 
#SBATCH --mem-per-cpu=15G
#SBATCH --mail-user=emily.beasley@umontreal.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 julia/1.9.1 
julia --project RunTheSim.jl 
