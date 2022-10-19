#!/bin/bash

# provide the maximum number of requested time for computation (in minutes)
#SBATCH -t 480

# number of requested compute nodes
#SBATCH --nodes=1

# load the necessary R module(s)
module load R/4.0.2

# change into the directory with the R code for the sim study
cd /R_scripts/

# run the R scripts with the sim study code
# Supply process identifier as argument
Rscript main_analyses.R ${SLURM_ARRAY_TASK_ID}

