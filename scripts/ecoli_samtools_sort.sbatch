#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=ecoli_sort_%J_stdout.txt
#SBATCH --error=ecoli_sort_%J_stderr.txt
#SBATCH --job-name=ecoli_sort
# 

bash /home/mbtoomey/BIOL7263_Genomics/scripts/BWA/ecoli_samtools_sort.sh