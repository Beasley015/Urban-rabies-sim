using DataFrames
using CSV
using DataAPI

# Write csv of parameters
serovals = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
barrier_vals = [0,1,2,3,4]
imm_sero = [0, 0.2, 0.4, 0.6, 0.8]
imm_disease = [0, 0.025, 0.05, 0.075, 0.1]


param_frame = DataAPI.allcombinations(DataFrame, "seros"=>serovals, "barrier"=>barrier_vals, "imm_sero"=>imm_sero, 
                                        "imm_disease"=>imm_disease)

CSV.write("params.csv", param_frame)

# Write the SLURM file 
job_file = """
#! /bin/bash
#SBATCH --array=1-$(size(param_frame, 1))
#SBATCH --time=02:30:00
#SBATCH --cpus-per-task=1 
#SBATCH --account=ctb-tpoisot
#SBATCH --output=$(joinpath("slurm", "%x-%a.out")) 
#SBATCH --mem-per-cpu=15G
#SBATCH --mail-user=emily.beasley@umontreal.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 julia/1.9.1 
julia --project RunTheSim.jl 
"""
write("rabies_ABM.sh", job_file)

# Use this line for full array size:
#$(size(param_frame, 1))