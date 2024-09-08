#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 8G
#SBATCH --output=pseud_fastqc_%J_stdout.txt
#SBATCH --error=pseud_fastqc_%J_stderr.txt
#SBATCH --job-name=ecoli_fastqc
# 

bash /home/mbtoomey/BIOL7263_Genomics/scripts/fastqc/pseud_fastqc.sh