#!/bin/bash


## Genome Reference files for GRCh37 to be downloaded prior building reference libraries.

BASE=${PWD}
GENOMES=${BASE}/Genomes

mkdir ${GENOMES}

## reference genome GRCh37
wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz -O ${GENOMES}/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

## reference annotation GRCh37
wget https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz -O ${GENOMES}/Homo_sapiens.GRCh37.87.gtf.gz
