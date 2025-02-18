Firstly, we need some to access some common software and download some data to help us on our way. For the software we will use a program called '[conda](https://docs.conda.io/en/latest/)':mag:, which will allow us to easily install lots of common bioinformatics software in a special 'environment' (think of it like a box :package:) without the need for admin access or other complications. For the data, we will download this from some public repositories.

Firstly, we need to enter or create a directory called "BIOL7263_Genomics". All further commands will be run within this  directory.



### Software
You will access software through the 'environment' :package: that I have set up for the class. This allows us to keep all the software in one place for easy access and repeatability (e.g. you may wish to run different versions of software for other analyses, you can do that in other environments). To access our class environment you will need to initate [mamba](https://mamba.readthedocs.io/en/latest/) and activate the environment by running the following at the OSCER command line


```bash
module load Mamba

mamba init

```
* Then log out, and log back into OSCER. This sign-out is required since Mamba modified your .bashrc file to make sure Mamba is set up properly right after the next time you log in.
* Once you logged back into OSCER, you should see your terminal beginning with `(base)`. This means conda `base` environment has been initialized successfully.

Now activate the class environment: 
```bash
mamba activate /home/mbtoomey/.conda/envs/BIOL7263_Genomics
```
Now you should see your terminal beginning with `(BIOL7263_Genomics)`. This means conda `BIOL7263_Genomics` environment has been initialized successfully. :warning: You must activate every time you log in. 


### Data
We will need to retrieve several sets of data for our adventure, this is similar to how you may collate data for your own analyses.
 1) Sequence Data
  * Either directly from a Sequencing Service or from a public access database.
 2) Reference Data
  * If you are lucky to have a reference genome...
 3) Databases
  * PFam-A

We will be working with the bacterial species *Escherichia coli* as it is a relatively small genome (which makes it easy for the timings of this tutorial), but the techniques you will learn here can be applied to any smaller or larger, and/or Eukaryotic genomes too!

#### Sequencing Data
Back at your home institute you will likely retrieve your data from either the institute's sequencing service or a private outside provider. However, there is also a wealth :moneybag: of sequenced genomic data stored in publically accesible repositories such as NCBI's [SRA](https://www.ncbi.nlm.nih.gov/sra) or EMBL-EBI's [ENA](https://www.ebi.ac.uk/ena). These portals are also where you will be required to deposit your sequencing efforts during publication.

For this adventure we will be downloading and processing raw sequencing data. Please note that some sequencing services may provide trimmed or quality assessed reads as part of their standard service, however it is up to you whether you want to use that data directly or process the raw data yourself. Always ask: are their methods directly suited to your analysis needs?

The raw data that we will use for the *E. coli* genome is available from [NCBI](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=ERR2789854) or [EMBL-EBI](https://www.ebi.ac.uk/ena/data/view/ERR2789854) with the accession ERR2789854 (they are also archived at the [DDBJ-DRA](https://www.ddbj.nig.ac.jp/dra/index-e.html) too). This is the same data but it is mirrored between the sites, however each site has a different way of accessing the data. We will focus on NCBI and EMBL-EBI for now.

With NCBI you need to use a tool called '[fastq-dump](https://ncbi.github.io/sra-tools/fastq-dump.html)':mag:, which given an accession and several other options will download the 'fastq' data files - it can be notoriously tricky and difficult at times and has some issues with downloading PacBio data. You can give it a try below if you wish, however the EMBL-EBI downloads will be much faster for this tutorial, so we strongly suggest you start there.

At EMBL-EBI they provide direct links to the 'fastq' files that were submitted to the archive ("Submitted files (FTP)"), and so you can use tools such as 'wget' or 'curl' to retrieve them.

NB - These commands may take a little bit of time to complete depending on your connection (NCBI: ~15  minutes; EMBL-EBI: ~2 minutes), so you might want to skip ahead to the next chapter for some light reading about sequencing technologies and file formats whilst you wait... don't forget to come back soon!

### Raw sequence data

We will download the data to our scratch folders. Let first make a directory called `BIOL7263_Genomics` and then make some sub directories within that:

```bash
cd scatch/[your user name]

mkdir BIOL7263_Genomics

cd BIOL7263_Genomics

mkdir -p sequencing_data/ecoli 

cd sequencing_data/ecoli

#download raw seq. 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/SRR857279/SRR857279_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/SRR857279/SRR857279_2.fastq.gz

# make the files read-only so we don't destroy our data accidentally
chmod 444 *.gz

# Now do the same for Chapter 5's Pseudomonas data

# Let's go back up to the main BIOL7263_Genomics folder
cd ../.. 

#then make a new folder
mkdir pseudomonas_gm41 

cd pseudomonas_gm41

# get the Illumina Data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR491/SRR491287/SRR491287_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR491/SRR491287/SRR491287_2.fastq.gz

# get the PacBio data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/006/SRR1042836/SRR1042836_subreads.fastq.gz

# make the files read-only so we don't destroy our data accidentally
chmod 444 *.gz
```

#### Reference Data
We will access the reference data from the National Center for Biotechnology Information (NCBI), check out the links below by Ctrl or Cmd clicking the link to open in a new tab:

[*Escherichia coli* str. K-12 substr. MG1655](https://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=161521)

There is a lot of information on this page, but the main pieces of information we are interested in are; the genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format, and the gene annotations in [GFF](https://en.wikipedia.org/wiki/General_feature_format) format. Can you see where these are? :eyes:

We will now download the data, as we are working with the command line we have already copied the links to the data below for you :slightly_smiling_face:. Using the '[wget](https://www.gnu.org/software/wget/)':mag: command we can download files directly from the web to our local dicrectory. The files are '*gzipped*', this means they are compressed to save space, it also allows us to make sure the data has not been corrupted during the transfer. We will also need to *unzip* them with the program '[gunzip](https://linux.die.net/man/1/gunzip)':mag:.

```bash
# Let's go back up to the main BIOL7263_Genomics folder
cd ..

# Create a directory to store our references
mkdir reference_sequences

cd reference_sequences

# Download the E. coli reference genome in FASTA and GFF formats
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz

# Make an ecoli directory within the references folder, and move the files there, and then unzip them
mkdir ecoli 

mv *.gz ecoli

gunzip ecoli/*.gz

# Change write permissions, so that we can't edit them by accident
chmod -R 444 ecoli/*.fna
chmod -R 444 ecoli/*.gff
```

#### Databases
We will need to get the PFam-A database of Hidden Markov Models (HMMS) and an Active Site database from the [Pfam](https://pfam.xfam.org/) website. They are located in an ftp directory. Use the commands below. Make sure you are in the "BIOL7263_Genomics" directory.

```bash
# Let's go back up to the main BIOL7263_Genomics folder
cd ../..

# create a directory and a sub-directory and move there
mkdir -p db/pfam && cd db/pfam

# Download the HMMs and .dat files needed for Pfam-A
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz

# Uncompress the files
gunzip *.gz
```

When the downloads are finished chack in with Prof. Toomey to double check your file structure. Then you may continue on to the adventure by clicking the title below.

# [Adventure Time!](https://github.com/mbtoomey/genomics_adventure/blob/release/chapter_2/task_1.md)
