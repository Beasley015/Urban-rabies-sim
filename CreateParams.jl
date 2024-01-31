using DataFrames
using CSV
using DataAPI

# Write csv of parameters
serovals = [0, 0.2, 0.4, 0.6, 0.8]
imm_type = ["wave", "propagule"]
barrier_vals = [0,1,2,3,4]
imm_sero = [0, 0.2, 0.4, 0.6, 0.8]
imm_disease = [0, 0.025, 0.05, 0.075, 0.1]


param_frame = DataAPI.allcombinations(DataFrame, "seros"=>serovals, "imm_type"=>imm_type, "barrier"=>barrier_vals, "imm_sero"=>imm_sero, 
                                        "imm_disease"=>imm_disease)

CSV.write("params.csv", param_frame)

# Write the SLURM file 
job_file = """
#! /bin/bash
#SBATCH --array=1-5
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1 
#SBATCH --output=$(joinpath("slurm", "%x-%a.out")) 

module load StdEnv/2020 julia/1.9.1 
julia --project RunTheSim.jl 
"""
write("rabies_ABM.sh", job_file)

# Use this line for full array size:
#$(size(param_frame, 1))