# Task 2 - QC the Data!

It is always important to check and understand the quality of the data you are working with.

Referring back to your earlier scripts (Chapter 2), write and execute .sbatch and .sh to run fastqc on the sequencing files in the pseudomonas_gm41 folder. My files are located here: `/scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/`

DOwnolad and open the files SRR1042836_subreads.fastq, SRR491287_1.fastq & SRR491287_2.fastq and check the reports generated.

![FastQC](https://github.com/mbtoomey/genomics_adventure/blob/release/images/chapter_5_task_2_image_1.png)

Note that the quality of the PacBio reads (SRR1042836_subreads.fastq) is much lower than the Illumina reads with a greater than 1 chance in 10 of there being a mistake for most reads.

![FastQC](https://github.com/mbtoomey/genomics_adventure/blob/release/images/chapter_5_task_2_image_2.png)

However, importantly, the length of the PacBio reads is much longer. Yay! :smile:.

# Go to [Task 3](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_5/task_3.md)
