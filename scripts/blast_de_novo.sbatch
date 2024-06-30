#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=blast_de_novo_%J_stdout.txt
#SBATCH --error=blast_de_novo_%J_stderr.txt
#SBATCH --job-name=blast_de_novo
# 

bash /home/mbtoomey/BIOL7263_Genomics/scripts/assembly/blast_de_novo.sh