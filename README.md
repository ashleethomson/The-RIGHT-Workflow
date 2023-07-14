# The RIGHT Workflow

This is an ALL Fusion calling pipeline that uses `Nextflow` (DSL2).

## Installation

Clone the repository to your system:

```bash
$ git clone --recurse-submodules https://github.com/ashleethomson/The-RIGHT-Workflow.git
```

## Download Reference genomes
To reproduce the results from "Reproducible bioinformatics analysis workflows for detecting Immunoglobulin Heavy Chain (IGH) gene fusions in B-cell acute lymphoblastic leukaemia patients", please run the `./DownloadGenome.sh` script in `nf-RIGHT-Workflow` to genome files used.

```bash
./DownloadGenome.sh
```

## Build Reference Libraries


Run the `BuildReferenceLibraries.sh` script to build the libraries for Arriba, STAR-Fusion, and FusionCatcher.

```bash
./BuildReferenceLibraries.sh
```

## Usage

### 1. Create a Sample-sheet

Create a _CSV_ file with the following structure:

```text
path,group,sample,filename,R1,R2
<path>,GROUP,TEST-0001,TEST-0001-XT,TEST-0001-XT_1.fastq.gz,TEST-0001-XT_2.fastq.gz
```

The columns are as follows:

- **path** = Path to parent directory of all ALL Fastq files
- **group** = Age group that sample belongs to e.g. _AYAII_
- **sample** = Sample identifier e.g. _AYAII_0001_
- **filename** = Fastq file basename e.g. _AYAII-0001-DIA1-PB_
- **R1/R2** = Fastq filename (e.g `AYAII-0001-DIA1-PB_R1.fastq.gz`)

The `group` and `sample` values are sub directories of `path`. The pipeline gets samples
by building the following `path` internally:

```text
|--------------- path --------------| |-group-|  |--sample--|
/cancer/storage/raw_fastq/ALL/RNAseq/   AYAII/   AYAII_0001/   <fastq reads at this level>
```

### 2. Run the Pipeline

Use the following command to run the pipeline:

```shell
nextflow run \
    main.nf \
    -profile slurm,conda \
    -N email@email.com \
    --outdir ./outdir \
    --samplesheet <path>/test.csv \
    --email email@gmail.com \
    --partition sahmri_prod_hpc
```
