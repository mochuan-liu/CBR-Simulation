#!/bin/bash

### Change to your own slurm setting
#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -n 25
#SBATCH -t 9-

### Change to your own R module available
module add r/4.1.3

Rscript simulation_setting_I_linear.R --seed 25587 --N 400 --method un
