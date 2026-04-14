#! /bin/bash
#SBATCH --array=1-450
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=1 
#SBATCH --account=def-tpoisot
#SBATCH --output=/home/ebeasley/slurm/%x-%a.out 
#SBATCH --mem-per-cpu=150G
#SBATCH --mail-user=ebeasley@bu.edu
#SBATCH --mail-type=ALL

module load StdEnv/2023 julia/1.9.3
julia --project RunTheSim.jl $SLURM_ARRAY_TASK_ID
