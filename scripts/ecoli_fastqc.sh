mkdir -p /scratch/mbtoomey/BIOL7263_Genomics/fastqc_output

fastqc /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/read_1.fastq.gz -o /scratch/mbtoomey/BIOL7263_Genomics/fastqc_output/
fastqc /scratch/mbtoomey/BIOL7263_Genomics/sequencing_data/ecoli/read_2.fastq.gz -o /scratch/mbtoomey/BIOL7263_Genomics/fastqc_output/