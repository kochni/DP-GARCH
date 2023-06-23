#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=15:00:00

module load gcc/8.2.0 r/4.2.2
Rscript R-REAL.R