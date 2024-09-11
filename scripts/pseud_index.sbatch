#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=pseud_index_%J_stdout.txt
#SBATCH --error=pseud_index_%J_stderr.txt
#SBATCH --job-name=pseud_index
# 

bash /home/mbtoomey/BIOL7263_Genomics/scripts/pseudomonas/pseud_index.sh