# Task 15 Locating Genes that are Missing Compared to the Reference
We can use a command from the [BEDTools](http://bedtools.readthedocs.org/en/latest/):mag: package to identify annotated genes which are not covered by reads across their full length.

For example, lets try the following:
```bash
bedtools coverage \
-a /scratch/mbtoomey/BIOL7263_Genomics/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.gff \
-b /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/ecoli_mapped_namesort_fixmate_sort_markdup.bam > /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/gene_coverage.txt
```

Here are the files I created: 
* [ecoli_cover.sh](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_cover.sh)
* [ecoli_cover.sbatch](https://github.com/mbtoomey/genomics_adventure/blob/release/scripts/ecoli_cover.sbatch)

This may take 10-15 mins once the job statrts. Take a look and you will see that the output contains one row per annotated gene, whereby the 13th (and final) column contains the proportion of the gene that is covered by reads from our sequencing. 1.00 means the gene is 100% covered and 0.00 means there is no coverage. 

Therefore, if we 'sort' our data by the 13th column we can see which genes are missing:
```bash
sort -t $'\t' -g -k 13 /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/gene_coverage.txt | less -S
```

It may be easier to view the data ommiting the annotation column, thus:

```bash
sort -t $'\t' -g -k 13 /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/mapping_to_reference/gene_coverage.txt | cut -f1-8,10-13 | less -S
```

For extra work, can you find these missing genes in IGV?

# Congratulations! :tada:
That concludes the first part of the adveture. You have successfully, Quality Controlled, filtered, mapped and analysed a whole bacterial genome! Well done! :trophy: :trophy: Time for a well deserved break! :coffee: :cookie:

In the next installment we will be looking at how to extract and assemble unmapped reads. This will enable us to look at material which may be present in the strain of interest but not in the reference sequence. See you soon! :wave:

# Now Start [Chapter 3](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_3/task_3.md)
