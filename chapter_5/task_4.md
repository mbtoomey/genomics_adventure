# Task 4 - Assembly Time!

Firstly let's construct an assembly using only the short-read Illumina data in a new directory!

```bash
cd /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/

mkdir Assembly 

cd Assembly
```

The next command, like most of the assemblies in this adventure (and in real life) take some time to complete.

```bash
spades.py --phred-offset 33 --threads 24 --careful -o /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/illumina_only \
-1 /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_trimmed_reads_val_1.fq.gz \
-2 /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_trimmed_reads_val_2.fq.gz
```
Create the .sh and .sbatch files locally and then upload and execute on OSCER. If your job fails, try increasing `mem` to `48G` and partition to `64gb_24core` in the sbatch. 

Here is how I setup my job: 
* [pseud_short_assembly.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/pseud_short_assembly.sh)
* [pseud_short_assembly.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/pseud_short_assembly.sbatch)

Let's take a look at the assembly. You know the score! Let's run QUAST!!
```bash
quast.py --output-dir /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/illumina_only/quast /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/illumina_only/contigs.fasta
```
Create the .sh and .sbatch files locally and then upload and execute on OSCER. 

Let's review the `quastreport.txt`

```
cat /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/illumina_only/quast/report.txt
```

Here is my result. You may get slightly different results here, as SPAdes uses a random seed.

```
Assembly                    contigs
# contigs (>= 0 bp)         528
# contigs (>= 1000 bp)      117
# contigs (>= 5000 bp)      100
# contigs (>= 10000 bp)     89
# contigs (>= 25000 bp)     74
# contigs (>= 50000 bp)     51
Total length (>= 0 bp)      6692019
Total length (>= 1000 bp)   6618751
Total length (>= 5000 bp)   6580437
Total length (>= 10000 bp)  6497588
Total length (>= 25000 bp)  6248274
Total length (>= 50000 bp)  5443414
# contigs                   127
Largest contig              248605
Total length                6626613
GC (%)                      59.00
N50                         102161
N75                         64662
L50                         22
L75                         43
# N's per 100 kbp           0.00
```

# Go to [Task 5](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_5/task_5.md)

