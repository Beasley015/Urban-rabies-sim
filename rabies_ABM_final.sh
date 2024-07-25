#! /bin/bash
#SBATCH --array=108
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=1 
#SBATCH --account=ctb-tpoisot
#SBATCH --output=slurm\%x-%a.out 
#SBATCH --mem-per-cpu=150G
#SBATCH --mail-user=emily.beasley@umontreal.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 julia/1.9.1 
julia --project RunTheSim.jl 
