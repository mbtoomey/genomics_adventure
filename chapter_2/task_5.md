# Task 5 - Generating an Index
Before we can start aligning reads to a reference genome, the genome sequence needs to be 'indexed', this means sorting the genome into easily searched chunks (similar to a book's index).

Now let's change to our *E. coli* reference sequence directory, do you remember where it is? Can you get there in one command? Also let's remind ourselves what is actually in the directory.
```bash
cd /scratch/mbtoomey/BIOL7263_Genomics/reference_sequences/ecoli

ls -lath
```

![directory listing](https://github.com/mbtoomey/genomics_adventure/blob/release/images/ref_folder.png)

In this directory we have 2 files:
 1. GCF_000005845.2_ASM584v2_genomic.fna - which is a FASTA file that contains the reference genome sequence.
 2. GCF_000005845.2_ASM584v2_genomic.gff - which is a file that contains the annotation for this genome.
 
You can have a look at them with 'more' or 'head' if you like to see the format they are in.

## BWA

Today we are using "BWA", but you may have read that is it deprecated/superceded by "minimap2" and there is a nice [blog](https://lh3.github.io/2018/04/02/minimap2-and-the-future-of-bwa)🔍 about it here. However, BWA works better with short (~100bp) reads and that is the data we have today. If you are feeling advanced you might want to run "minimap2" instead through the rest of the adventure, but you're heading off into Mirkwood without help...Watch out for orcs! 

Now type 'bwa' in your terminal and see what happens. Hopefully, you should see something similar to this:

![bwa help](https://github.com/mbtoomey/genomics_adventure/blob/release/images/bwa_help.png)

BWA is actually a suite of programs which perform many different functions. We are only going to use two of these functions during our adventure, 'bwa index' and 'bwa mem'. Explore the outputs of these two commands now.

By default 'bwa index' will use the '*IS*' algorithm to produce the index. This works well for most genomes, but for very large ones (e.g. vertebrates) you may need to use '*bwtsw*' algorithm. For bacterial genomes the default algorithm will work fine.

### Create a Reference Index

We will submit this as a job so let's first make a folder for our scripts in our home folder: 

```bash
mkdir /home/mbtoomey/BIOL7263_Genomics/scripts/BWA/
```

Then make a .sh file with the following: 

```bash
bwa index /scratch/mbtoomey/BIOL7263_Genomics/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna -p /scratch/mbtoomey/BIOL7263_Genomics/reference_sequences/ecoli/
```

Here the -p option sets the output folder, so our index will be written to the same location as our refrence sequence. Now upload the .sh and the .sbatch files the scripts/BWA folder and submit the job. 

Here are my scripts: 
* [ecoli_trim.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_trim.sh)
* [ecoli_trim.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_trim.sbatch)

Take a look at the output of this command in your terminal (e.g. 'ls -lath'). You will notice that the 'bwa index' program has created a set of new files (e.g .sa, .amb, etc). We don't need to know what these are right now, just that these are the index files that 'bwa mem' will need.

![BWA index output](https://github.com/mbtoomey/genomics_adventure/blob/release/images/bwa__help_index_output.png)

# Go to [Task 6](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_2/task_6.md)
