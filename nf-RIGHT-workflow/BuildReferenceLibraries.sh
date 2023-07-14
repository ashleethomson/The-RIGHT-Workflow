#!/bin/bash

## Reference files for each algorithm need to be built prior to running the pipeline

BASE=${PWD}
GENOMES=${BASE}/Genomes
REF_LIB=${BASE}/ReferencesLibraries
STARINDEX=${REF_LIB}/STAR_Index_GRCh37
REF_ARRIBA=${REF_LIB}/Arriba
REF_FC=${REF_LIB}/FusionCatcher
REF_STARFUSION=${REF_LIB}/STAR_Fusion

mkdir ${REF_LIB}
mkdir ${STARINDEX}
mkdir ${REF_ARRIBA}
mkdir ${REF_FC}
mkdir ${REF_STARFUSION}



## STAR Index

STAR --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir ${STARINDEX} \
    --genomeFastaFiles ${GENOMES}/GRCh37.dna.primary_assembly.fa.gz \
    --sjdbGTFfile ${GENOMES}/Homo_sapiens.GRCh37.87.gtf.gz \
    --sjdbOverhang 74


## Arriba References

wget https://github.com/suhrig/arriba/releases/download/v2.1.0/arriba_v2.1.0.tar.gz -O ${REF_ARRIBA}/arriba_v2.1.0.tar.gz
    tar -xzvf ${REF_ARRIBA}/arriba_v2.1.0.tar.gz
    rm ${REF_ARRIBA}/arriba_v2.1.0.tar.gz
    mv ${REF_ARRIBA}/arriba_v2.1.0/database/* ${REF_ARRIBA}/
    rm -r ${REF_ARRIBA}/arriba_v2.1.0


## FusionCatcher

wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.aa -O ${REF_FC}/human_v102.tar.gz.aa
wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ab -O ${REF_FC}/human_v102.tar.gz.ab
wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ac -O ${REF_FC}/human_v102.tar.gz.ac
wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ad -O ${REF_FC}/human_v102.tar.gz.ad
cat ${REF_FC}/human_v102.tar.gz.* > ${REF_FC}/human_v102.tar.gz
rm -f ${REF_FC}/human_v102.tar.gz.*
tar xvf ${REF_FC}/human_v102.tar.gz
mv ${REF_FC}/human_v102.tar.gz/* ${REF_FC}/



## STAR-Fusion

wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz --no-check-certificate -O ${REF_STARFUSION}/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz
tar xvf ${REF_STARFUSION}/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz
rm ${REF_STARFUSION}/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz
mv ${REF_STARFUSION}/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir ${REF_STARFUSION}/
rm -f ${REF_STARFUSION}/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play

