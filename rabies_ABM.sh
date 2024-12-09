#! /bin/bash
#SBATCH --array=106,109,117,118,119,120,137,15,16,17,18,182,19,191,198,20,201,206,21,215,216,22,269,274,305,320,340,341,342,384,396,398,399,400,401,404,414,415,416,83,84,85
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1 
#SBATCH --account=ctb-tpoisot
#SBATCH --output=slurm\%x-%a.out 
#SBATCH --mem-per-cpu=100G
#SBATCH --mail-user=emily.beasley@umontreal.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 julia/1.9.1 
julia --project RunTheSim.jl 
