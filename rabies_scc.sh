#!/bin/bash -l

# Set SCC project
#$ -P dietzelab

# Request buyin nodes so Mike doesn't kill you
#$ -l buyin

# Specify array job with tasks
#$ -t 1-450

# Specify hard time limit of hours for the job (if you have a short runtime the SCC gives you priority)
#$ -l h_rt=20:00:00

# Assign cores and cores per node
#$ -pe omp 28

# Send an email when the job finishes or if it is aborted 
#$ -m n
#
#
# Below is what would get passed to the command line - you can test just these lines to make sure they work

cd /projectnb/dietzelab/ebeasley/UrbanRabies

module load julia/1.9.1

julia /projectnb/dietzelab/ebeasley/UrbanRabies/RunTheSim.jl $SGE_TASK_ID
