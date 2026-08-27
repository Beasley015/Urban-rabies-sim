#! /bin/bash
#SBATCH --array=1-27
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1 
#SBATCH --account=def-tpoisot
#SBATCH --output=slurm\%x-%a.out 
#SBATCH --mem-per-cpu=100G
#SBATCH --mail-user=ebeasley@bu.edu
#SBATCH --mail-type=ALL

module load StdEnv/2023 julia/1.9.3 
julia --project ParamSensitivity.jl $SLURM_ARRAY_TASK_ID
