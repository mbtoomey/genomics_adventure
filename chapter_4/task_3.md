# Map Reads Back to Assembly
Here we will use BWA again to index the contigs.fasta file and remap the reads. This is almost identical to the procedure we followed during the alignment section, the only difference is that instead of aligning to the reference genome, we are aligning to our newly created reference.

Make sure that you are in the directory "/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly". Let's create a new directory to keep our work separate and organized/

```bash
mkdir mapping_to_assembly

cd mapping_to_assembly
```

We now need to index and map our files, this is similar to the steps from Chapter 2... Indeed we will use the QC reads we created in Task 2! Here are the steps one-by-one, but let's try and combine them into a single script for submission. 

```bash
# Index the contigs
bwa index /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/contig.fasta

# align QC reads to contigs and output SAM file
bwa mem -t 2 /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/contigs.fasta \
/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_1.fq.gz \
/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_1.fq.gz \
> /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/contigs_mapped.sam
```

Once complete we can convert the SAM file to a BAM file:
```bash
samtools view -bS /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/contigs_mapped.sam >  /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/contigs_mapped.bam
```

And then we can sort the BAM file:
```bash
samtools sort -o /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/contigs_mapped_sorted.bam \
/scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/contigs_mapped.bam
```

Once completed, we can index the BAM file:

```
samtools index /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/contigs_mapped_sorted.bam
```

<details>
  <summary>Advanced: All in one command</summary>
```bash
  bwa index /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/contigs.fasta
bwa mem -t 2 /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/contigs.fasta /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_1.fq.gz /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_2.fq.gz > /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/contigs_mapped.sam 
samtools sort -o /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/contigs_mapped_sorted.bam /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/contigs_mapped.sam
samtools index /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/contigs_mapped_sorted.bam
```

Here my scripts for this job:

* [align_de_novo.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/align_de_novo.sh)
* [align_de_novo.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/align_de_novo.sbatch)

</details>

We can then (at last!) obtain some basic summary statistics using the samtools flagstat command:
```
cd /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/

samtools flagstat contigs_mapped_sorted.bam

5903641 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
1581 + 0 supplementary
0 + 0 duplicates
5903421 + 0 mapped (100.00% : N/A)
5902060 + 0 paired in sequencing
2951030 + 0 read1
2951030 + 0 read2
0 + 0 properly paired (0.00% : N/A)
5901840 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
394 + 0 with mate mapped to a different chr
2 + 0 with mate mapped to a different chr (mapQ>=5)
```

We can see that very few of the reads do not map back to the contigs. Importantly 99% of the reads are properly paired, which gives us some indication that there are not too many mis-assemblies.

We can run 'qualimap' to get some more detailed information (and some images too), it'll take a couple of minutes:

```bash
qualimap bamqc -outdir /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/bamqc -bam /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/contigs_mapped_sorted.bam
```

<details>
  <summary>Here my scripts for this job:</summary>
* [qmap_de_novo.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/qmap_de_novo.sh)
* [qmap_de_novo.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/qmap_de_novo.sbatch)

<details>


Download the `bamqc` folder and open the qualimapReport.html file in your browser. 

Go to either the "Chromosome stats" section and you will see that the larger of our contigs have a mean coverage of around 160 - which is what we would expect from our original alignment.

If you notice very carefully :wink:, there is one contig which has a size of 46899 - this is very very close to the size (46850) of the main contig we found in the unmapped reads assembly - another good indication that it is a separate sequence (remember we suspected it was a plasmid) and not integrated into a chromosome. We can double check this with a quick blast search...

```bash
blastn -subject /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/contigs.fasta \
-query /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli//unmapped_assembly/spades_assembly/contigs.fasta \
-outfmt 6 -out /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/mapping_to_assembly/check_plasmid.blastn
```

Opening the 'check_plasmid.blastn' we can see the top hit as:
```
NODE_1_length_46850_cov_69.441152       NODE_33_length_46899_cov_70.549934      100.000 46850   0       0       1       46850   50      46899   0.0     86516
```

This shows us that this contig almost exactly matches the the unmapped assembly, strongly supporting that this is a plasmid sequence and not integrated into the chromosomes.

# Go to [Task 4](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_4/task_4.md)
