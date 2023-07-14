#!/bin/bash


######## Must update job-name and output with date MANUALLY ############

#SBATCH --job-name=TheRIGHTworkflow
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=/path/to/log/directory/TheRIGHTworkflow.log

## Resources allocation request parameters
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=40000       # In Mb, this will request 80 Gb of memory
#SBATCH --time=36:00:00          # Max run time in hh:mm:s


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