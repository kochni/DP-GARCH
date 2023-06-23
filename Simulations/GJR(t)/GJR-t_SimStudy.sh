#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=41
#SBATCH --mem-per-cpu=3G
#SBATCH --time=2:00:00

module load gcc/8.2.0 r/4.2.2
Rscript GJR-t_SimStudy.R