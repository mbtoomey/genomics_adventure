# convert to bam
samtools view -bS /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseud_pacbio.sam > /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio.bam

# sort
samtools sort -o /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio_sorted.bam /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio.bam

# index
samtools index /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio_sorted.bam

# stats
samtools flagstat /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio_sorted.bam > /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio_sorted.stats