#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 8G
#SBATCH --output=ecoli_trim_%J_stdout.txt
#SBATCH --error=ecoli_trim_%J_stderr.txt
#SBATCH --job-name=ecoli_trim
# 

bash /home/mbtoomey/BIOL7263_Genomics/scripts/trim_galore/ecoli_trim.sh