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
cat /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/assembly/quast/quastreport.txt

cat /scratch/mbtoomey/BIOL7263_Genomics/pseudomonas_gm41/hybrid/quast/quastreport.txt
```

Here is my result. You may get slightly different results here, as SPAdes uses a random seed.
```
Assembly                    illumina_only_contigs  hybrid_contigs
# contigs (>= 0 bp)         612                    265           
# contigs (>= 1000 bp)      323                    94            
# contigs (>= 5000 bp)      253                    84            
# contigs (>= 10000 bp)     205                    79            
# contigs (>= 25000 bp)     88                     62            
# contigs (>= 50000 bp)     29                     41            
Total length (>= 0 bp)      6642052                6684539       
Total length (>= 1000 bp)   6589395                6658678       
Total length (>= 5000 bp)   6395291                6632315       
Total length (>= 10000 bp)  6038360                6593847       
Total length (>= 25000 bp)  4113769                6281030       
Total length (>= 50000 bp)  2097787                5518221       
# contigs                   339                    96            
Largest contig              187340                 671305        
Total length                6601889                6660346       
GC (%)                      58.98                  58.97         
N50                         32206                  137540        
N75                         19064                  62570         
L50                         60                     15            
L75                         127                    32
```

Hooray! It seems that using the longer reads has improved the completeness of our assembly - reducing the number of contigs and increasing the N50.

# Go to [Task 6](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_5/task_6.md)
