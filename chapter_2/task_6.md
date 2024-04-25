# Task 6 - Mapping Reads to the Indexed Reference Sequence
Now we can begin to align 'read_1' and 'read_2' to the reference genome. First of all change back into raw reads directory where we made the trimmed and QC reads, and then create a subdirectory to contain our mapping results.
```bash
cd /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/

mkdir mapping_to_reference

cd mapping_to_reference
```

## Mapping
Have a look at the output from typing 'bwa mem'. You should see something like this:

![bwa mem](https://github.com/mbtoomey/genomics_adventure/blob/release/images/bwa_mem_help.png)

The basic format of the command show as:

>Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]

From this we can see that we need to provide BWA with a set of FASTQ files containing the raw reads (denoted by
'\<in.fq>' meaning required and '[in2.fq]' as optional) and a reference file (unhelpfully this is listed as '\<idxbase>') and any other options we wish to change. 
  
Our reference sequences are in the file:
>/scratch/mbtoomey/BIOL7263_Genomics/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna

Our filtered reads are in the files:
>/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_1.fq.gz

>/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_2.fq.gz

So, in order to align our paired reads using multi-threading, and output to a file called ecoli_mapped.sam. We will set the output '-o'  to ***ecoli_mapped.sam*** and the number of threads/processes the program should use ('-t'). The multithreading option allows for parallel processing which should speed up the run. However, to take advantage of the we will need to specify it in our .sh file (below) and also in the 'n-tasks' options of the .sbatch file 

```bash
bwa mem -t 4 \
/scratch/mbtoomey/BIOL7263_Genomics/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic \
/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_1.fq.gz \
~/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_2.fq.gz
-o /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/ecoli_mapped.sam
```

Now create the .sh files with the above commands updated for your file system and a .sbatch file with the 'n-task = 4' option, upload to the scripts/BWA folder and submit the job.

Here are the files I created: 
* [ecoli_bwa_mem.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_bwa_mem.sh)
* [ecoli_bwa_mem.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_bwa_mem.sbatch)

This process will take a couple of minutes to complete. You can monitor the progress by reviewing the stdout.txt files, These will look something like this: 

![bwa mapping](https://github.com/mbtoomey/genomics_adventure/blob/release/images/bwa_mem_run.png)

Congratulations, you have just performed your first mapping of reads to a reference genome!! :trophy: But that only get's us so far... There's always more to do.

## SAM File Manipulation - Basic
Once the alignment is complete, list the directory contents and check that the alignment file 'ecoli_mapped.sam' is present. The raw alignment is stored in what is called 'SAM' format (Simple AlignMent format). It is a plain text file, and so you can view it using the 'less', 'head', or 'tail' commands, for example. It is best to not open the whole file in a text editor, as you will likely run out of memory, they can get very big! In our case it is about 2.2GB.

Each alignment line has 11 mandatory fields for essential alignment information including mapping position, and a variable number of optional fields for flexible or aligner specific information. For further details as to what each field means see the PDF [here](http://samtools.sourceforge.net/SAM1.pdf) :mag:. If you are on the Workshop on Genomics course then Mike Zody will have explained these files to you in gory detail! We will also explore this more in the next task.


# Go to [Task 7](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_2/task_7.md)
