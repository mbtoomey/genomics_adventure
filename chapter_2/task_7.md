# Task 7 - SAM File Manipulation - A Bit Less Basic...
Before we can visualize the alignment in a more meaningful way, it is beneficial to convert the SAM file into a BAM file (Binary AlignMent format) for several reasons. SAM files are great for human readability, but they are not so good for fast computational access. The binary format will allow speedy access to the information stored within it, and it also reduces the file size on the disk.

To do this we will use another suite of programs called '[samtools](http://www.htslib.org/)' :mag:. Go ahead and type the command below, and have a look at the options:
```bash
samtools view
```

## Converting SAM to BAM
We can see that we need to provide 'samtools view' with a reference genome as a FASTA formatted file '-T', the '-b' and '-S' options to indicate that the output should be in BAM format, and that the input is in SAM format, all along with the alignment input file itself - e.g.
```bash
samtools view -b -S \
-T /scratch/mbtoomey/BIOL7263_Genomics/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped.sam \
-o /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped.bam
```

Now create the .sh file with the above commands updated for your file system and a .sbatch file, upload to the scripts/BWA folder and submit the job.

Here are the files I created: 
* [ecoli_samtools_view.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_samtools_view.sh)
* [ecoli_samtools_view.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_samtools_view.sbatch)

NB - Again we are using the '-o' option, but a simple redirect is also possible. You may also see the options contracted as '-bS' in other tutorials, and it is perfectly acceptable to do that, but in the beginning I like to make it clear we are using two separate options so we know exactly what is going on.

This should take a few minutes. Note, that for larger datasets you may wish to set multiple threads as well with the --threads option. Remember it you will also need to specify the number of threads in the .sbatch file with the "ntasks" option. 

It's always good to check that your files have processed correctly, if something goes wrong it's better to
catch it now rather than ten commands later. Have you noticed that the 'BAM' file (712Mb) is much smaller than the 'SAM' file (2.2GB) - this is to be expected as the binary format is more efficient at storing information - think of zip files as a similarish concept.

## Sorting a BAM File
Once the conversion is complete we will need to sort the BAM file, this is so that the reads are stored in the order they appear along the 'chromosomes', or 'scaffolds', or 'contigs' in your reference. We can do this using the 'samtools sort' command. There are a couple of different ways to sort a BAM file, and it is left to the user to decide which is best. For the next task we will use the default options and it will take about 2 minutes to complete.

```bash
samtools sort /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped.bam \
-o /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped_sorted.bam
```

## A Note on Efficiency
In the previous set of tasks (4, 5 & 6) we have aligned the trimmed reads to the reference genome, converted the SAM to a BAM and then sorted the resulting BAM file. For clarity, and the purposes of learning we have done this in individual steps. However, in 'real-life' :tm:, it is much faster and easier to do these steps in one single command using Unix pipes! There is no need to do the next steps, of course, as we have already done them above but the commands are given below in order for you to compare.

For example, these are the three commands we completed previously:
```bash
bwa mem -t 4 \
/scratch/mbtoomey/BIOL7263_Genomics/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic \
/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_1.fq.gz \
~/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_2.fq.gz \
-o /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped.sam

samtools view -b -S \
-T /scratch/mbtoomey/BIOL7263_Genomics/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped.sam \
-o /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped.bam

samtools sort /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped.bam \
-o /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped_sorted.bam
```

but we can also write it as one command, like this:
```bash
bwa mem -t 4 \
/scratch/mbtoomey/BIOL7263_Genomics/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic \
/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_1.fq.gz \
/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_2.fq.gz \
-o /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped.sam
| samtools sort -O bam -o /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped_sorted_onecommand.bam
```
Notice how in the second example we have used one less command too. This is because the output of 'bwa mem' is as a SAM file, and therefore it can be piped directly into the 'samtools sort' which is able to output in BAM format.

We can put these together into a job submission and save time because we will not need to wait for multiple jobs to queue. 

Here are the files I created: 
* [ecoli_map_sort.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_map_sort.sh)
* [ecoli_map_sort.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_map_sort.sbatch)

# Go to [Task 8](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_2/task_8.md)
