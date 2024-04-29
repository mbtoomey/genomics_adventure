#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 8G
#SBATCH --output=ecoli_qualmap_%J_stdout.txt
#SBATCH --error=ecoli_qualmap_%J_stderr.txt
#SBATCH --job-name=ecoli_qualmap
# 

bash /home/mbtoomey/BIOL7263_Genomics/scripts/BWA/ecoli_qualmap.sh