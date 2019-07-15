# Download base image ubuntu 16.04
FROM ubuntu:16.04
# Update Ubuntu Software repository
RUN apt-get update
# Copy BRCA-analyzer
COPY brca_analyzer/ /brca_analyzer/
# Unzip reference files of BRCA-analyzer
CMD gunzip /brca_analyzer/ref/*.gz
# Install python3
RUN apt-get -y install python3 python3-pip default-jdk python3-biopython less vim
# Install python3-packages
RUN pip3 install argparse numpy python-docx xlrd xlsxwriter matplotlib
# Install BRCA-analyzer
RUN cd /brca_analyzer/ && python3 install.py -bwa /bwa/ -sam /samtools/ -picard /picard/dist/ -gatk /GATK/ -snpeff /snpEff/ -annovar /annovar/
# Install ANNOVAR
COPY annovar/ /annovar/
# Install cpanminus
RUN apt-get -y install cpanminus make wget
# Install perl module
RUN cpanm Pod::Usage
# Install samtools
COPY samtools-1.2/ /samtools/
# Install bwa
COPY bwa-0.7.15/ /bwa/
# Install picard
COPY picard-2.0.1/ /picard/
# Install GATK
COPY GenomeAnalysisTK-3.7/ /GATK/
# Install snpEff
COPY snpEff/ /snpEff/
# Unzip archives for ANNOVAR
CMD gunzip /annovar/humandb/*.gz
