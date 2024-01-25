using DataFrames
using CSV
using DataAPI

serovals = [0, 0.2, 0.4, 0.6, 0.8]
imm_type = ["wave", "propagule"]
barrier_vals = [0,1,2,3,4]
imm_sero = [0, 0.2, 0.4, 0.6, 0.8]
imm_disease = [0, 0.025, 0.05, 0.075, 0.1]


param_frame = DataAPI.allcombinations(DataFrame, "seros"=>serovals, "imm_type"=>imm_type, "barrier"=>barrier_vals, "imm_sero"=>imm_sero, 
                                        "imm_disease"=>imm_disease)

CSV.write("params.csv", param_frame)