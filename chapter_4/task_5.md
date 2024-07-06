# Task 5 - Annotation of *de novo* Assembled Contigs
We will now annotate the contigs using BLAST and Pfam as with the unmapped contigs.

As before, we’ll 'call' open reading frames within the *de novo* assembly. We will use codon table 9 (transl_table=11) which defines the bacterial [codon usage table](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)🔍. We will also restrict the ORFs to just those sequences longer than 300 nucleotides (i.e. 100 amino acids). We will generate two files a fasta "orfs.fa" with the predicted amino acid sequences and a .bed file with the genome locations for mapping in IGV. 

```bash
orfipy /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/contigs.fasta --pep orfs.fa --min 10 --max 10000 --procs 4 --min 300 --table 9 --outdir /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly

orfipy /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/contigs.fasta --bed /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/orfs.bed --table 9 --min 300 --procs 4
```
<details>
  <summary>Here my scripts for this job:</summary>
  
* [orfipy.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/orfipy.sh)
* [orfipy.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/orfipy.sbatch)

</details>

As with the unmapped reads we will search the open reading frames against the Pfam HMM database of protein families. Later on we will be able to use these results to identify Pfam domains which are unique to a particular strain.

This will take around 1 hour on the OSCER, so you may want to submit this job and return to this later 

```bash
pfam_scan.pl -fasta /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/contigs.orf.fasta \
-dir  /scratch/mbtoomey/BIOL7263_Genomics/db/pfam/ \
-outfile /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/contigs.orf.pfam -cpu 2 -as
```
<details>
  <summary>Here my scripts for this job:</summary>
  
* [pfam.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/pfam.sh)
* [pfam.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/pfam.sbatch)

</details>

Another simple way to annotate the genes is to blast the orfs.fa file against known proteins. You could blast all possible proteins in the "nr" database, but this would take a large amount of memory and time. We can speed up the search by creating a smaller database of proteins from an existing [E. coli genome assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/). 

```bash
cd /scratch/mbtoomey/BIOL7263_Genomics/db

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_translated_cds.faa.gz

gunzip GCF_000005845.2_ASM584v2_translated_cds.faa.gz

makeblastdb -in GCF_000005845.2_ASM584v2_translated_cds.faa -dbtype prot -parse_seqids -out Ec_prot
```

We can then use blastp to find the top hit in this smaller database for each putative reading frame in the genome: 

```bash
blastp -query /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/orfs.fa -db /scratch/mbtoomey/BIOL7263_Genomics/db/Ec_prot -outfmt "6 qseqid sseqid pident stitle" -max_target_seqs 1 | sort -u > /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/assembly/orf_hit.txt
```
Were we have limited the output to a single hit `-max_target_seqs 1` and output the percentage identity (pident) and peptide name (stitle).  [This page](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6) has a useful summary of the output options for command line blast. 

<details>
  <summary>Here my scripts for this job:</summary>
  
* [orf_blast.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/orf_blast.sh)
* [orf_blast.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/orf_blast.sbatch)

</details>

We can then merge the blast output with the orf bed file to add names to the bed file. To do this I simply downloaded the `orfs.bed` and `orf_hit.txt` and then joined and rearranged with a few tidyverse commands in R: 

* [annotate_bed.R](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/bed_file.R)

Now you can load the resulting bed file (`orfs_gene.bed`) to the contigs.fasta in IGV. 

![igv](https://github.com/mbtoomey/genomics_adventure/blob/release/images/chapter_4_task_5_image_1.png)

# [Chapter 5](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_5/task_1.md) awaits..!
