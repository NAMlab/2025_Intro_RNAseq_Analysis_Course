# Here we get the data using a script from ENA Browser, checks the qualitiy of our fastq files and performs an alignment to the human cDNA transcriptome reference with Kallisto.
# This script checks the qualitiy of our fastq files and performs an alignment to the human cDNA transcriptome reference with Kallisto.
# To run this 'shell script' you will need to open your terminal and navigate to the directory where this script resides on your computer.
# This should be the same directory where you fastq files and reference fasta file are found.
# Change permissions on your computer so that you can run a shell script by typing: 'chmod +x script.sh' (without the quotes) at the terminal prompt 
# Then type './script.sh' (without the quotes) at the prompt.  
# This will begin the process of running each line of code in the shell script.

#Generate the directories
mkdir data
mkdir qc
mkdir genomes
mkdir mapped


SRR11855627
SRR11855628
SRR11855629
SRR11855630
SRR11855631
SRR11855632


# Get the data (using ENA https://www.ebi.ac.uk/ena/browser/home) 
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/027/SRR11855627/SRR11855627_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/027/SRR11855627/SRR11855627_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/028/SRR11855628/SRR11855628_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/028/SRR11855628/SRR11855628_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/029/SRR11855629/SRR11855629_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/029/SRR11855629/SRR11855629_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/030/SRR11855630/SRR11855630_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/030/SRR11855630/SRR11855630_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/031/SRR11855631/SRR11855631_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/031/SRR11855631/SRR11855631_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/032/SRR11855632/SRR11855632_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/032/SRR11855632/SRR11855632_2.fastq.gz



# Get the data (using SRA & sra-tools https://www.ncbi.nlm.nih.gov/sra/)
fasterq-dump --split-files -O data SRR11855627 && gzip data/SRR11855627_*.fastq
fasterq-dump --split-files -O data SRR11855628 && gzip data/SRR11855628_*.fastq
fasterq-dump --split-files -O data SRR11855629 && gzip data/SRR11855629_*.fastq
fasterq-dump --split-files -O data SRR11855630 && gzip data/SRR11855630_*.fastq
fasterq-dump --split-files -O data SRR11855631 && gzip data/SRR11855631_*.fastq
fasterq-dump --split-files -O data SRR11855632 && gzip data/SRR11855632_*.fastq



# To direct the output to a specific subdirectory use: wget -nc -P output_subdirectory input_file 


# first use fastqc to check the quality of our fastq files:
fastqc data/*.gz -t 4 -o qc

# trimm your reads if needed after the qc check: http://www.usadellab.org/cms/?page=trimmomatic
# trimmomatic PE -threads 8 input1 input2 paired_output1 unpaired_output1 paired_output2 unpaired_output2 \
# ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 \
# HEADCROP:5 \
# SLIDINGWINDOW:4:15 \
# MINLEN:30

trimmomatic PE -threads 8 data/SRR11855627_1.fastq.gz data/SRR11855627_2.fastq.gz data/SRR11855627_1_trimmed_paired.fastq.gz data/SRR11855627_1_trimmed_unpaired.fastq.gz data/SRR11855627_2_trimmed_paired.fastq.gz data/SRR11855627_2_trimmed_unpaired.fastq.gz ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic PE -threads 8 data/SRR11855628_1.fastq.gz data/SRR11855628_2.fastq.gz data/SRR11855628_1_trimmed_paired.fastq.gz data/SRR11855628_1_trimmed_unpaired.fastq.gz data/SRR11855628_2_trimmed_paired.fastq.gz data/SRR11855628_2_trimmed_unpaired.fastq.gz ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic PE -threads 8 data/SRR11855629_1.fastq.gz data/SRR11855629_2.fastq.gz data/SRR11855629_1_trimmed_paired.fastq.gz data/SRR11855629_1_trimmed_unpaired.fastq.gz data/SRR11855629_2_trimmed_paired.fastq.gz data/SRR11855629_2_trimmed_unpaired.fastq.gz ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic PE -threads 8 data/SRR11855630_1.fastq.gz data/SRR11855630_2.fastq.gz data/SRR11855630_1_trimmed_paired.fastq.gz data/SRR11855630_1_trimmed_unpaired.fastq.gz data/SRR11855630_2_trimmed_paired.fastq.gz data/SRR11855630_2_trimmed_unpaired.fastq.gz ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic PE -threads 8 data/SRR11855631_1.fastq.gz data/SRR11855631_2.fastq.gz data/SRR11855631_1_trimmed_paired.fastq.gz data/SRR11855631_1_trimmed_unpaired.fastq.gz data/SRR11855631_2_trimmed_paired.fastq.gz data/SRR11855631_2_trimmed_unpaired.fastq.gz ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic PE -threads 8 data/SRR11855632_1.fastq.gz data/SRR11855632_2.fastq.gz data/SRR11855632_1_trimmed_paired.fastq.gz data/SRR11855632_1_trimmed_unpaired.fastq.gz data/SRR11855632_2_trimmed_paired.fastq.gz data/SRR11855632_2_trimmed_unpaired.fastq.gz ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30


# next, we want to build an index from our reference fasta file 
# I get my reference mammalian transcriptome files from here: https://useast.ensembl.org/info/data/ftp/index.html
kallisto index -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index genomes/Solanum_lycopersicum.SL3.0.cdna.all.fa.gz

# first use fastqc to check the quality of the trimmed fastq files:
fastqc data/*trimmed.fastq.gz -t 4 -o qc

# now map reads to the indexed reference transcriptome (cdna) 


kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR11855627 -t 8 data/SRR11855627_1_trimmed_paired.fastq.gz data/SRR11855627_2_trimmed_paired.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR11855628 -t 8 data/SRR11855628_1_trimmed_paired.fastq.gz data/SRR11855628_2_trimmed_paired.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR11855629 -t 8 data/SRR11855629_1_trimmed_paired.fastq.gz data/SRR11855629_2_trimmed_paired.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR11855630 -t 8 data/SRR11855630_1_trimmed_paired.fastq.gz data/SRR11855630_2_trimmed_paired.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR11855631 -t 8 data/SRR11855631_1_trimmed_paired.fastq.gz data/SRR11855631_2_trimmed_paired.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR11855632 -t 8 data/SRR11855632_1_trimmed_paired.fastq.gz data/SRR11855632_2_trimmed_paired.fastq.gz

# summarize fastqc and kallisto mapping results in a single summary html using MultiQC
multiqc -d . 
 
echo "Finished"