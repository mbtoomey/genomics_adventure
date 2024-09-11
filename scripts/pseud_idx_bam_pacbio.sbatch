#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=pseud_idxbam_PB_%J_stdout.txt
#SBATCH --error=pseud_idxbam_PB_%J_stderr.txt
#SBATCH --job-name=pseud_idxbam_PB
# 

bash /home/mbtoomey/BIOL7263_Genomics/scripts/pseudomonas/pseud_idx_bam_pacbio.sh