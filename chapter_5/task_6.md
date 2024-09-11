# Task 6 - Align Reads Back to Reference
Let's realign our original reads back to the assembly and see what we have - refer to previous notes if you are unsure of the steps.

First create an index with BWA from the hybrid assembly.
```bash
#make a new folder for the mapping output
mkdir /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/

#make the index and write it to the new folder
bwa index -p /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/hybrid_assembly \
/scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/hybrid/contigs.fasta 
```

First map the Illumina reads, and then follow the standard protocols we learned previously...You can follow along, or make your own commands...

```bash
# map
bwa mem -t 6 /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/hybrid_assembly \
/scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_trimmed_reads_val_1.fq.gz \
/scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_trimmed_reads_val_2.fq.gz \
-o /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseud_illumina.sam

# convert to bam
samtools view -bS /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseud_illumina.sam > /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina.bam

# sort
samtools sort -o /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina_sorted.bam /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina.bam

# index
samtools index /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina_sorted.bam

# stats
samtools flagstat /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina_sorted.bam > /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina_sorted.stats
```

We can now map our PacBio data to our assembly too!. For this we will use another tool called "minimap2" which is better suited to mapping PacBio data than BWA. Read more on this [here](https://lh3.github.io/2018/04/02/minimap2-and-the-future-of-bwa) :mag:
It too creates SAM files, and then we can follow the same procedure with the output as we do for Illumina data.
```bash
minimap2 -x map-pb -t 6 -a /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/hybrid/contigs.fasta \
/scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR1042836_subreads.fastq.gz  \
-o /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseud_pacbio.sam

# convert to bam
samtools view -bS /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseud_pacbio.sam > /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio.bam

# sort
samtools sort -o /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio_sorted.bam /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio.bam

# index
samtools index /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio_sorted.bam

# stats
samtools flagstat /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio_sorted.bam > /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_pacbio_sorted.stats
```

Have a look at both flagstat files and give them a comparison.

# Go to [Task 7](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_5/task_7.md)
