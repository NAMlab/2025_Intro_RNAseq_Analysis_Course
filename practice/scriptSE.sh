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

# Get the data (using ENA https://www.ebi.ac.uk/ena/browser/home) 
wget -nc -P data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/003/SRR2006793/SRR2006793.fastq.gz
wget -nc -P data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/007/SRR2006797/SRR2006797.fastq.gz
wget -nc -P data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/002/SRR2006792/SRR2006792.fastq.gz
wget -nc -P data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/004/SRR2006794/SRR2006794.fastq.gz
wget -nc -P data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/001/SRR2006791/SRR2006791.fastq.gz
wget -nc -P data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/005/SRR2006795/SRR2006795.fastq.gz
wget -nc -P data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/006/SRR2006796/SRR2006796.fastq.gz
wget -nc -P data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/008/SRR2006798/SRR2006798.fastq.gz


# Get the data (using SRA & sra-tools https://www.ncbi.nlm.nih.gov/sra/)
fasterq-dump --split-files -O data SRR2006793 && gzip data/SRR2006793*.fastq
fasterq-dump --split-files -O data SRR2006797 && gzip data/SRR2006797*.fastq
fasterq-dump --split-files -O data SRR2006792 && gzip data/SRR2006792*.fastq
fasterq-dump --split-files -O data SRR2006794 && gzip data/SRR2006794*.fastq
fasterq-dump --split-files -O data SRR2006791 && gzip data/SRR2006791*.fastq
fasterq-dump --split-files -O data SRR2006795 && gzip data/SRR2006795*.fastq
fasterq-dump --split-files -O data SRR2006796 && gzip data/SRR2006796*.fastq
fasterq-dump --split-files -O data SRR2006798 && gzip data/SRR2006798*.fastq


# To direct the output to a specific subdirectory use: wget -nc -P output_subdirectory input_file 


# first use fastqc to check the quality of our fastq files:
fastqc data/*.gz -t 4 -o qc

# trimm your reads if needed after the qc check: http://www.usadellab.org/cms/?page=trimmomatic
# trimmomatic SE -threads 8 input output \
# ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 \
# HEADCROP:5 \
# SLIDINGWINDOW:4:15 \
# MINLEN:30

trimmomatic SE -threads 8 data/SRR2006793.fastq.gz data/SRR2006793_trimmed.fastq.gz ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic SE -threads 8 data/SRR2006797.fastq.gz data/SRR2006797_trimmed.fastq.gz ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic SE -threads 8 data/SRR2006792.fastq.gz data/SRR2006792_trimmed.fastq.gz ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic SE -threads 8 data/SRR2006794.fastq.gz data/SRR2006794_trimmed.fastq.gz ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic SE -threads 8 data/SRR2006791.fastq.gz data/SRR2006791_trimmed.fastq.gz ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic SE -threads 8 data/SRR2006795.fastq.gz data/SRR2006795_trimmed.fastq.gz ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic SE -threads 8 data/SRR2006796.fastq.gz data/SRR2006796_trimmed.fastq.gz ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30
trimmomatic SE -threads 8 data/SRR2006798.fastq.gz data/SRR2006798_trimmed.fastq.gz ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30

 
# next, we want to build an index from our reference fasta file 
# I get my reference mammalian transcriptome files from here: https://useast.ensembl.org/info/data/ftp/index.html
kallisto index -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index genomes/Solanum_lycopersicum.SL3.0.cdna.all.fa.gz

# first use fastqc to check the quality of the trimmed fastq files:
fastqc data/*trimmed.fastq.gz -t 4 -o qc
 
# now map reads to the indexed reference transcriptome (cdna)

kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR2006793 -t 4 --single -l 50 -s 30 data/SRR2006793_trimmed.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR2006797 -t 4 --single -l 50 -s 30 data/SRR2006797_trimmed.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR2006792 -t 4 --single -l 50 -s 30 data/SRR2006792_trimmed.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR2006794 -t 4 --single -l 50 -s 30 data/SRR2006794_trimmed.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR2006791 -t 4 --single -l 50 -s 30 data/SRR2006791_trimmed.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR2006795 -t 4 --single -l 50 -s 30 data/SRR2006795_trimmed.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR2006796 -t 4 --single -l 50 -s 30 data/SRR2006796_trimmed.fastq.gz
kallisto quant -i genomes/Solanum_lycopersicum.SL3.0.cdna.all.index -o mapped/SRR2006798 -t 4 --single -l 20 -s 30 data/SRR2006798_trimmed.fastq.gz
 

 
 
# summarize fastqc and kallisto mapping results in a single summary html using MultiQC
multiqc -d . 
 
echo "Finished"
