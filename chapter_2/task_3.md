# Task 3 - Random Subsampling and Digital Normalisation
Sometimes you may find that you have sequenced too much data! :open_mouth: However, this is not as bad a place to be in as it first seems, e.g it may put limits on your infrastructure; for example, your computer/server/HPC may not be able to keep all the data in memory to assemble it. Or you may want to do some quick analyses to sanity check your data without doing a full analysis.

There are several methods of reducing your data, here we will discuss two of them:

### 1. Random Subsampling
As easy as the name suggests, we take a random subsample of the original dataset, e.g. 10% of the data and then we can use that data to perform an assembly, although we should really use all the original data to do any mapping. Why do you think this might be the case?

For this task we will use the program ['seqtk'](https://github.com/lh3/seqtk) :mag: which is an excellent little toolkit to do FASTA/Q processing.

Check it out:
```bash
seqtk

seqtk sample
```

Now let's try randomly subsampling our *E. coli* 'read_1' dataset to 10% of the data. We need to express 10% as a fraction for 'seqtk', so 0.1, and supply that as an option to the 'seqtk sample' program. Let's do it twice, to see if we truly do get a random subsampling of the data...
```bash
seqtk sample trimmed_reads_val_1.fq.gz 0.1 > trimmed_reads_val_1_subsample_one.fq

seqtk sample trimmed_reads_val_1.fq.gz 0.1 > trimmed_reads_val_1_subsample_two.fq

head trimmed_reads_val_1_subsample*
```

What do the headers look like, do they look like a random selection to you? Ha ha! :stuck_out_tongue_closed_eyes: I played another trick on you! Don't worry though, this is a common pitfall when you start out learning to use bioinformatic software.

So, what is going on? Well 'seqtk sample', by default, sets a special 'seed' value of '11' - it is sneakily there in the help section - and so in this case both files happen to be the same random subsampled selection! :dizzy_face:

This may seem somewhat counter-intertuitive at first - how can something be random if you can repeat it exactly - however it is immensely useful, as you will come to see. The program uses a special trick (called [Reservoir sampling](https://en.wikipedia.org/wiki/Reservoir_sampling) :mag: - you don't need to know this, just to be aware) that will take the same 'random' sample when given a starting 'seed' - in this case '11' - allowing it to start the random selection process at the same place each time. Change the 'seed', and you will get a different subsampled selection! Neat. :sunglasses:

Using the '-s' option we can set the 'seed' starting position of the random subsampling. So, if we change the value to '1234' and '5678' we should have two different random samples. For this next bit you can choose any numbers you like! :metal:
```bash
seqtk sample -s 1234 trimmed_reads_val_1.fq.gz 0.1 > trimmed_reads_val_1_subsample_three.fq

seqtk sample -s 5678 trimmed_reads_val_1.fq.gz 0.1 > trimmed_reads_val_1_subsample_four.fq

head trimmed_reads_val_1_subsample_three.fq

head trimmed_reads_val_1_subsample_four.fq
```

How do they look now? Much more random right!? :thumbsup: Of course, if you are doing this in groups (or on multiple machines) the random samples should be the same give the same seed.

Let's tidy up our testing files now, and then repeat the process one last time - we will need this next set of files later on in our adventure. Please use the 'seed' value that is indicated, or your analyses will look very different to the examples shown here.
```bash
rm *.fq

seqtk sample -s 628 trimmed_reads_val_1.fq.gz 0.5 > trimmed_reads_val_1_subsampled.fq

seqtk sample -s 628 trimmed_reads_val_2.fq.gz 0.5 > trimmed_reads_val_2_subsampled.fq
```

There is just one final thing we need to do, although 'seqtk' can read 'gzipped' files, it does not produce them. So we should 'gzip' them now. You can either use the program 'gzip' or it's speedier cousin 'pigz'.
```bash
pigz *.fq
```

Now your sequence directory should look something like this: 

![subsamples](https://github.com/mbtoomey/genomics_adventure/blob/release/images/subsamples.png)

One thing to be careful of when randomly subsampling data like this is that we may, by accident, subsample regions in our genome that are more, or less, sequenced than others - especially if our coverage is not even across the sequencing libraries. This can lead to issues with downstream analysis, especially if for example, the libraries come from Multiple Displacement Amplification (e.g. Single Cells). So we should always be careful with this approach. One method to minimise this potential source of bias in our data is called Digital Normalisation.

### 2. Digital Normalisation

In all seriousness, I forgot to write this section. Mea culpa. So as it's 2024 and we are all using ChatGPT anyway, I gave it the prompt "write me a few short paragraphs about digital normalisation of fastq illumina reads".

Let's have a look at what it came up with.

>Certainly! Let's delve into the concept of **digital normalization** for **FASTQ Illumina reads**.

>### Digital Normalization of FASTQ Illumina Reads

>**Digital normalization** is a computational technique used to reduce redundancy and improve the efficiency of downstream analysis for high-throughput sequencing data, particularly **FASTQ files** generated by Illumina sequencers. Here are the key points:

>1. **FASTQ Format**:
   - **FASTQ** is a widely adopted text-based file format that stores both raw sequence data and quality scores.
   - Illumina sequencing systems produce data in FASTQ format, which includes the actual DNA sequences (reads) and their corresponding quality scores.
   - These quality scores reflect the confidence in the base calls at each position in the read.

>2. **Redundancy and Data Volume**:
   - When dealing with large-scale sequencing projects, the sheer volume of data can be overwhelming.
   - Many reads are similar or identical due to PCR amplification, sequencing errors, and other factors.
   - Redundant reads contribute to increased storage requirements, computational time, and memory usage during downstream analysis.

>3. **Digital Normalization Process**:
   - Digital normalization aims to retain a representative subset of reads while discarding redundant ones.
   - It involves identifying and removing duplicate or nearly identical reads based on sequence similarity.
   - By doing so, the data volume is reduced without compromising the overall information content.

>4. **Benefits**:
   - **Storage Efficiency**: Digital normalization significantly reduces the size of FASTQ files, making them more manageable for storage and transfer.
   - **Computational Speed**: Smaller datasets lead to faster processing times during alignment, variant calling, and other analyses.
   - **Quality Preservation**: The technique retains high-quality reads while discarding low-quality or repetitive ones.

>5. **Downstream Analysis**:
   - After digital normalization, the remaining reads can be used for various downstream analyses, such as **alignment** against a reference genome.
   - Alignment tools map the short reads to specific genomic locations, allowing researchers to study genetic variations, gene expression, and more.

>In summary, digital normalization is a valuable preprocessing step that optimizes the use of sequencing data, streamlines subsequent analyses, and facilitates meaningful biological insights. Researchers can choose from various tools and pipelines to perform this normalization, tailoring the approach to their specific research goals¹[1] ²[2].

>Source: Conversation with Bing, 4/25/2024
(1) Sequence File Formats | FASTQ & BCL formats for Illumina sequencing. https://support.illumina.com/content/illumina-marketing/en/informatics/sequencing-data-analysis/sequence-file-formats.html.
(2) Alignment – NGS Analysis. https://learn.gencore.bio.nyu.edu/alignment/.
(3) Simulated Illumina BRCA1 reads in FASTQ format - figshare. https://f1000.figshare.com/articles/dataset/Simulated_Illumina_BRCA1_reads_in_FASTQ_format/92338.
(4) NGS : reads quality control - INRA. https://urgi.versailles.inra.fr/download/tutorial/NGS-reads-quality-control/Galaxy-reads-quality-control.pdf.
(5) NGS Data Analysis for Illumina Platform—Overview and Workflow. https://www.thermofisher.com/us/en/home/life-science/cloning/cloning-learning-center/invitrogen-school-of-molecular-biology/next-generation-sequencing/ngs-data-analysis-illumina.html.


Not too bad! We won't use this technique today but a good tool is [BBNorm](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbnorm-guide/) 🔍 from the BBTools package.

## Contaminant Checking
A number of tools are available which also enable to you to quickly search through your reads and assign them to particular taxa or taxonomic groups. These can serve as a quick check to make sure your samples or libraries are not contaminated with DNA from other sources. If you are performing a de-novo assembly, for example, and have DNA sequences present from multiple organisms, you will risk poor results and chimeric contigs.

Some ‘contaminants’ may turn out to be inevitable by-products of sampling and DNA extraction, and this is often the case with algae, and/or other symbionts but some groups have made amazing discoveries such as the discovery of a third symbiont (which turned out to be a yeast) in lichen, see [here](http://science.sciencemag.org/content/353/6298/488.full) :mag:.

We won't cover this topic in our current adventure, but here is a list of some tools you can use to check the taxonomic classification of reads:
 * [Kraken](https://ccb.jhu.edu/software/kraken2/)
 * [Centrifuge](https://ccb.jhu.edu/software/centrifuge/)
 * [Blobology](https://blobtoolkit.genomehubs.org/)
 * Blast (in conjunction with sub-sampling your reads) and Krona to plot results
 * and many more!

# Go to [Task 4](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_2/task_4.md)
