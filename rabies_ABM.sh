#! /bin/bash
#SBATCH --array=1-5
#SBATCH --time=28:00:00
#SBATCH --cpus-per-task=1 
#SBATCH --account=ctb-tpoisot
#SBATCH --output=/home/ebeasley/slurm/%x-%a.out 
#SBATCH --mem-per-cpu=150G
#SBATCH --mail-user=ebeasley@bu.edu
#SBATCH --mail-type=ALL

module load StdEnv/2020 julia/1.9.1 
julia --project RunTheSim.jl 
