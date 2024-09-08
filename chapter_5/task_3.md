# Task 3 - Trim the Reads!

Trim the Illumina short-reads as your did in chapter 2. Here is my command

```bash
trim_galore --paired --gzip --cores 4 --length 100 \
/scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_1.fastq.gz \
/scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_2.fastq.gz \
--basename SRR491287_trimmed_reads -o /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/
```

Create the .sh and .sbatch files locally and then upload and execute on OSCER.

You can check the number of filtered reads using "grep –c" and the quality of trimmed reads with fastqc if you want.

For our next trick we want to keep the long reads from PacBio even though they are of lower quality. We are relying
on the assembler to use them appropriately...

# Go to [Task 4](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_5/task_4.md)
