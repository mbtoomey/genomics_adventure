# Chapter Two - QC, Alignment and Variant Calling
## Task 1a - Evaluating the Quality of Illumina Data
From your terminal (command line), navigate to the 'sequencing_data/ecoli' directory (you may be there already) and list the contents of the directory.
```bash
cd /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli

ls -lath
```
![Terminal showing listing of fastq.gz files](https://github.com/mbtoomey/genomics_adventure/blob/release/images/ecoli_folder_peek.png)

These are the sequencing reads from a porject sequencing the genome of Escherichia coli Bal225 with paired-end 2x250bp Illumina MiSeq that I downloaded from the [SRA](https://www.ncbi.nlm.nih.gov/sra/SRX282099[accn]). The filenames are a bit long and boring to continue to type out - so let's make them a little easier by making a symbolic link to them (this is a type of special file that acts as a pointer to the real file).
```bash
ln -s SRR857279_1.fastq.gz read_1.fastq.gz
ln -s SRR857279_2.fastq.gz read_2.fastq.gz
```
Now is we peek at the files again: 

```bash
ls -lath
```
![with symbolic links](https://github.com/mbtoomey/genomics_adventure/blob/release/images/ecoli_symbol_link.png)
The symbolic links appear as files in the folder and we can treat them link files in our further analyses

Now we have nice and easy file names to work with :smile:. You should always record or note down what your original filenames are, so that you can refer to the correct data in the future! Just like keeping notes in a lab book. :open_book:

Most programs that work with sequence data require that the 'read 1' and 'read 2' files have the reads in the same order. You can identify reads from the same pair because they will have the same header followed by either a "1" or a "2". We will now look at the raw reads to make sure they look 'correct'. To view the first few headers we can use the 'zcat' command (this is similar to 'cat' but works with zipped files), 'head' to see the top three lines, and then 'grep' to catch the header lines (in this case we know they start with "@SRR" - you may need to check your own files for a suitable header).
```bash
zcat read_1.fastq.gz | head | grep @SRR
zcat read_2.fastq.gz | head | grep @SRR
```

![head of fastq files](https://github.com/mbtoomey/genomics_adventure/blob/release/images/ecoli_grep.png)

The only difference in the headers for the two reads is the read number. Of course, this is no guarantee that all the headers in the file are consistent. To get some more confidence, lets repeat the above commands using 'tail' instead of 'head' to compare reads at the end of the files. This will take a little longer, as tail has to read all the way through to the end of the file!
```bash
zcat read_1.fastq.gz | tail | grep @SRR
zcat read_2.fastq.gz | tail | grep @SRR
```

You can also check that there is an identical number of reads in each file using 'zcat', 'grep' and 'wc -l'. This should take about 4 minutes for each.
```bash
zcat read_1.fastq.gz | grep @SRR | wc –l
zcat read_2.fastq.gz | grep @SRR | wc –l
```

Oops! :trollface: Did you just copy and paste that and receive and error saying "wc: –l: No such file or directory" :stuck_out_tongue_winking_eye:. That's okay! Remember, typing the commands so that you get used to using them, and so that you understand what the options do is a much better way of learning! Don't forget you can use 'Tab complete' to automatically complete filenames. Anyway, the answer you should have received is '4273258. Try again, with the correct command, and see what you get! :hugs: Hint: the "–" was wrong in the above command...

# Quality control

Now, let's run the 'fastqc' program to examine the quality of the raw reads. 

The 'fastqc' program performs a number of tests which determines whether a green tick (pass) :heavy_check_mark:, exclamation mark (warning) :exclamation:, or a red cross (fail) :x: is displayed. However, it is important to realise that fastqc has no knowledge of what your specific library is or should look like. All of its tests are based on a completely random library with 50% GC content. Therefore, if you have a sample which does not match these assumptions, it may 'fail' the library.

For example, if you have a high AT or high GC organism it may fail the 'per sequence GC content' test. If you have any barcodes or low complexity libraries (e.g. small RNA libraries, RAD-Seq, Amplicons) they may also fail some of the sequence complexity tests. The bottom line is that you need to be aware of what your library is and whether what 'fastqc' is reporting makes sense for that type of library. Know your data from study design to sequencing and beyond!

Since fastqc will require more computational resource, so will will submit this action as a job through using SLURM. First we will write a shell script (.sh file) will our commands and then an sbatch file to submit that script. 

First lets setup a folder in our home directory for the scripts we will make

```bash
mkdir -p /home/mbtoomey/BIOL7263_Genomics/scripts/fastqc
```

Now create a .sh file with the following commands and save it to ***/home/mbtoomey/BIOL7263_Genomics/scripts/fastqc***. This simplest way to do this is to create the file offline in an editor like [Notepad++](https://notepad-plus-plus.org/) and then upload to OSCER with [scp or winSCP](https://www.ou.edu/oscer/support/file_transfer). I created a file called ***ecoli_fastqc.sh*** that contains the following:

```bash
mkdir -p /scratch/mbtoomey/BIOL7263_Genomics/fastqc_output

fastqc /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/read_1.fastq.gz -o /scratch/mbtoomey/BIOL7263_Genomics/fastqc_output/
fastqc /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/read_2.fastq.gz -o /scratch/mbtoomey/BIOL7263_Genomics/fastqc_output/
```

Next we will need to create an sbatch file that the SLURM job manager will read to queue and start our job. I created a file called ***ecoli_fastqc.sbatch*** that contained the following:

```bash
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 8G
#SBATCH --output=ecoli_fastqc_%J_stdout.txt
#SBATCH --error=ecoli_fastqc_%J_stderr.txt
#SBATCH --job-name=ecoli_fastqc
# 

bash /home/mbtoomey/BIOL7263_Genomics/scripts/fastqc/ecoli_fastqc.sh
```
Here are the files I created:
* [ecoli_fastqc.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_fastqc.sh)
* [ecoli_fastqc.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_fastqc.sbatch)

Once the files are uploaded we simple submit the job with the following command:

```bash
sbatch home/mbtoomey/BIOL7263_Genomics/scripts/fastqc/ecoli_fastqc.sbatch
```

Oh no! I got an error: 

![endline error](https://github.com/mbtoomey/genomics_adventure/blob/release/images/line_break1.png)

This happened because I created the file in the windows program that uses a different encoding for the ends of the lines of text. In Notepad++ I can correct this with this option: 

![endline error correction](https://github.com/mbtoomey/genomics_adventure/blob/release/images/line_break2.png)

Now we if we uploaded the corrected scrpt and resubmit in will be accepted and our job is given an ID# and we can check on its progress with ***squeue -u mbtoomey*** 

![Job status](https://github.com/mbtoomey/genomics_adventure/blob/release/images/fastqc_submission.png)

Notice that the job is pending. This happens when there are many jobs on the cluster. We will have to wait for an indeterminate amount of time. Such is life when you are using a free resource :shrug:

...sometime later....

WHen the job runs it will produce stderror and stdoutput txt files in the folder where the job was submitted from. These provide detail of run progress and any error messages that usually apprear on screen when you are running in interactive mode. 

![error and out files](https://github.com/mbtoomey/genomics_adventure/blob/release/images/error_out.png)

We can uses ***cat*** to look at the out file: 

![out file](https://github.com/mbtoomey/genomics_adventure/blob/release/images/std_out.png)

Looks like everything ran just fine. Let's take a look at the fastqc output. This output is an html file so we will need to download the fastqc_output folder and examine locally. I used WInSCP to do this: 

![fastqc out](https://github.com/mbtoomey/genomics_adventure/blob/release/images/fastqc_out.png)

## FastQC output

![FastQC: Basic Stats](https://github.com/mbtoomey/genomics_adventure/blob/release/images/fcqc1.png)

In this case we have a number of errors and warnings which at first sight suggest that there has been a problem - but don't worry too much yet. Let's go through a few of them.

### Per base sequence quality
This is one of the most important metrics. If the quality scores are poor, either the wrong FASTQ encoding has been guessed by fastqc (see the Basic Statistics tab), or the data itself is poor quality. This view shows an overview of the range of quality values across all bases at each position in the FASTQ file.  Generally anything with a median quality score greater than Q20 is regarded as acceptable; anything above Q30 is regarded as 'good'. For more details, see the help documentation in fastqc.

![FastQC: Per base sequence quality](https://github.com/mbtoomey/genomics_adventure/blob/release/images/fcqc2.png)

In this case this check is red - and it is true that the quality drops off at the end of the reads. It is normal for read quality to get worse towards the end of the read. You can see that at ~100 bases the quality is still relatively good.

### Per base sequence Content
For a completely randomly generated library with a GC content of 50% one expects that at any given position within a read there will be a 25% chance of finding an A,C,T or G base. Here we can see that our library satisfies these criteria, although there appears to be some minor bias at the beginning of the read. This may be due to PCR duplicates during amplification or during library preparation. It is unlikely that one will ever see a perfectly uniform distribution.

![FastQC: Per base sequence Content](https://github.com/mbtoomey/genomics_adventure/blob/release/images/fcqc3.png)

### Sequence Duplication Levels
In a library that covers a whole genome uniformly most sequences will occur only once in the final set. A low level of duplication may indicate a very high level of coverage of the target sequence, but a high level of duplication is more likely to indicate some kind of enrichment bias (e.g. PCR over-amplification).
This module counts the degree of duplication for every sequence in the set and creates a plot showing the relative number of sequences with different degrees of duplication. 

![FastQC: Sequence Duplication](https://github.com/mbtoomey/genomics_adventure/blob/release/images/fcqc4.png)

### Overrepresented Sequences
This checks for sequences that occur more frequently than expected in your data. It also checks any sequences it finds against a small database of known sequences. A typical cause is that the original DNA was shorter than the length of the read - so the sequencing overruns the actual DNA and runs into the adaptors used to bind it to the flow cell.

![FastQC: Overrepresented Sequences](https://github.com/mbtoomey/genomics_adventure/blob/release/images/fcqc5.png)

### Adaptor Content
In this case it has found that a small number of reads (35000) that appear to contain a sequence used in the preparation for the library. Don't worry, as we can trim these in a later stage and is completely normal to find them in your data.

![FastQC: Adapter content](https://github.com/mbtoomey/genomics_adventure/blob/release/images/fcqc6.png)

### Other Reports
Have a look at them and at what the author of FastQC has to say [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/):mag:. Or check out their youtube tutorial video [here](https://www.youtube.com/watch?v=bz93ReOv87Y):mag:.

Remember the error and warning flags are the author's (albeit experienced) judgement of what typical data should look like. It is up to you to use some initiative and understand whether what you are seeing is typical for your dataset and how that might affect any analysis you are performing.

## Task 1b - Evaluating the Quality of Illumina Data Continued...
Do the same for 'read 2' as you did for 'read 1', (you can open a new file in the same window) and have a look at the various plots and metrics which are generated. How similar are they? Why might they differ?

Note that the number of reads reported in both files are identical. Overall, both 'read 1' and 'read 2' can be regarded as 'good enough' data-sets. 

For reference, [here](https://mbtoomey.github.io/genome_biology_FA24/fastqc_examples/36_F_S113_R1_001_fastqc.html) and [here](https://mbtoomey.github.io/genome_biology_FA24/fastqc_examples/41_G_S80_R1_001_fastqc.html) are a few "bad" dataset with extreme levels of duplication, adaptor content, and quality drop offs mid-sequence. These results suggest that something went wrong in library preparation :sob: 

## Now Go to [Task 2](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_2/task_2.md)
