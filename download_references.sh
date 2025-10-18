#!/bin/bash

mkdir genome_ref
cd genome_ref
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.p14.genome.fa.gz -O hg38.fa.gz
gunzip hg38.fa.gz
cd ..

mkdir epigenome_ref
cd epigenome_ref
mkdir K562
cd K562
wget https://www.encodeproject.org/files/ENCFF754EAC/@@download/ENCFF754EAC.bigWig -O ATAC.bigWig
wget https://www.encodeproject.org/files/ENCFF413AHU/@@download/ENCFF413AHU.bigWig -O DNase.bigWig
wget https://www.encodeproject.org/files/ENCFF849TDM/@@download/ENCFF849TDM.bigWig -O H3K27ac.bigWig
wget https://www.encodeproject.org/files/ENCFF253TOF/@@download/ENCFF253TOF.bigWig -O H3K4me3.bigWig
cd ..
cd ..
