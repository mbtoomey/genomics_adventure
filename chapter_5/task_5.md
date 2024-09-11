# Create a Hybrid Assembly
Now will execute the same command, but this time include the longer PacBio reads to see the effect it has on our assembly.

```bash
spades.py --phred-offset 33 --threads 24 --careful \
-o /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/hybrid \
--pacbio /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR1042836_subreads.fastq.gz \
-1 /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_trimmed_reads_val_1.fq.gz \
-2 /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_trimmed_reads_val_2.fq.gz
```
Here is how I setup my job: 
* [pseud_long_assembly.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/pseud_long_assembly.sh)
* [pseud_long_assembly.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/pseud_long_assembly.sbatch)

Well, you definitely know the score! Let's run QUAST!! 

```bash
quast.py --output-dir /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/hybrid/quast /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/hybrid/contigs.fasta
```

Now let's review and compare the `quastreport.txt` for the illumina and hybrid assemblies. 

```
cat /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/illumina_only/quast/report.txt

cat /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/hybrid/quast/report.txt
```

Here is my result. You may get slightly different results here, as SPAdes uses a random seed.
```
Assembly                    illumina_only_contigs  
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

Assembly                    hybrid_contigs
# contigs (>= 0 bp)         300
# contigs (>= 1000 bp)      32
# contigs (>= 5000 bp)      30
# contigs (>= 10000 bp)     28
# contigs (>= 25000 bp)     25
# contigs (>= 50000 bp)     20
Total length (>= 0 bp)      6725842
Total length (>= 1000 bp)   6677010
Total length (>= 5000 bp)   6672688
Total length (>= 10000 bp)  6656532
Total length (>= 25000 bp)  6595226
Total length (>= 50000 bp)  6437266
# contigs                   34
Largest contig              1292976
Total length                6678398
GC (%)                      58.98
N50                         489830
N75                         355188
L50                         4
L75                         8
# N's per 100 kbp           0.00
```

Hooray! It seems that using the longer reads has improved the completeness of our assembly - reducing the number of contigs and increasing the N50.

# Go to [Task 6](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_5/task_6.md)
