#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=ecoli_idxbam_%J_stdout.txt
#SBATCH --error=ecoli_idxbam_%J_stderr.txt
#SBATCH --job-name=ecoli_idxbam
# 

bash /home/mbtoomey/BIOL7263_Genomics/scripts/BWA/ecoli_idxbam.sh