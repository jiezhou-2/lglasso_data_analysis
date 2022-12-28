#!/bin/bash
#SBATCH --job-name=lglasso
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=real1.out
#SBATCH --error=real1.err
time R CMD BATCH network_real_data.R test1.out
