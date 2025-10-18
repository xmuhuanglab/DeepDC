#!/bin/bash

mkdir -p genome_ref
cd genome_ref
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.p14.genome.fa.gz -O hg38.fa.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
gunzip hg38.fa.gz
chmod -R 755 ../chopchop/bowtie/
../chopchop/bowtie/bowtie-build -f ./hg38.fa hg38 --threads 4
cd ..

mkdir -p epigenome_ref
cd epigenome_ref
mkdir -p K562
cd K562
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF754EAC/@@download/ENCFF754EAC.bigWig -O ATAC.bigWig
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF413AHU/@@download/ENCFF413AHU.bigWig -O DNase.bigWig
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF849TDM/@@download/ENCFF849TDM.bigWig -O H3K27ac.bigWig
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF253TOF/@@download/ENCFF253TOF.bigWig -O H3K4me3.bigWig
cd ..
cd ..
