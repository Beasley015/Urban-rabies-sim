using DataFrames
using CSV
using DataAPI

# Write csv of parameters
serovals = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
barrier_vals = [0,1,2,3,4]
imm_disease = [0, 0.5, 0.75, 0.1, 0.125]
imm_type = ["propagule","wave"]


param_frame = DataAPI.allcombinations(DataFrame, "seros"=>serovals, "barrier"=>barrier_vals, "imm_disease"=>imm_disease,
                                                "imm_type"=>imm_type)

CSV.write("params.csv", param_frame)

# Write the SLURM file 
job_file = """
#! /bin/bash
#SBATCH --array=1-$(size(param_frame, 1))
#SBATCH --time=09:00:00
#SBATCH --cpus-per-task=1 
#SBATCH --account=ctb-tpoisot
#SBATCH --output=$(joinpath("slurm", "%x-%a.out")) 
#SBATCH --mem-per-cpu=100G
#SBATCH --mail-user=emily.beasley@umontreal.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 julia/1.9.1 
julia --project RunTheSim.jl 
"""
write("rabies_ABM.sh", job_file)

# Use this line for full array size:
print(size(param_frame, 1))
