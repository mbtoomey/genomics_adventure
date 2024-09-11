#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem 16G
#SBATCH --output=pseud_bwa_%J_stdout.txt
#SBATCH --error=pseud_bwa_%J_stderr.txt
#SBATCH --job-name=pseud_bwa
# 

bash /home/mbtoomey/BIOL7263_Genomics/scripts/pseudomonas/pseud_bwa_mem.sh