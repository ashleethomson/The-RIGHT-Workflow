#!/bin/bash


###### BASE PARAMETERS STAY THE SAME BETWEEN RUNS #######
## Directories
BASE=/path/to/The-Right-Workflow
OUTDIR=${BASE}/outdir
SAMPLESHEET_DIR=${BASE}/samplesheets

## Variables
NF_MAIN=nf-RIGHT-workflow/main.nf
EMAIL=email@email.com

##### VARIABLE PARAMETERS CHANGE WITH EACH RUN #######
## DATE must be same format as samplesheet
DATE=<DATE>
QUEUE=<#>


nextflow run \
    ${BASE}/${NF_MAIN} \
    -profile slurm,conda \
    -N ${EMAIL} \
    --outdir ${OUTDIR}/${DATE} \
    --samplesheet ${SAMPLESHEET_DIR}/samplesheet_${DATE}.csv \
    --email ${EMAIL} \
    -qs ${QUEUE} \
    --genome /path/to/genome \
    --annotation /path/to/annotation \
    --workDir <str>
