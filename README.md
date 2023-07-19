# The RIGHT Workflow

This is an ALL Fusion calling pipeline that uses `Nextflow` (DSL2).

## Installation

Clone the repository to your system:

```bash
$ git clone https://github.com/ashleethomson/The-RIGHT-Workflow.git
```

## Reference libraries and genomes
The GRCh37 reference genome and annotation was used for "Reproducible bioinformatics analysis workflows for detecting _IGH_ gene fusions in B-cell acute lymphoblastic leukaemia patients". To replicate the results found in this paper, please use this version of the genome as well.   

Reference libraries for each algorithm will need to be downloaded from their respective sources for this pipeline, as seen in the `nextflow.config` file, and a STAR Index needs to be created. You will have to install STAR to create the Index. STAR version 2.7.9a was used in the manuscript.  

The STAR Index used was created using the following parameters:
```bash
STAR \
    --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir star-2.7.9a-75bp \
    --genomeFastaFiles GRCh37.fa \
    --sjdbGTFfile ref-transcripts.gtf \
    --sjdbOverhang 74
```

## Usage

### 1. Create a Sample-sheet

Create a _CSV_ file with the following structure:

```text
path,group,sample,filename,R1,R2
<path>,GROUP,TEST-0001,TEST-0001-XT,TEST-0001-XT_1.fastq.gz,TEST-0001-XT_2.fastq.gz
```

The columns are as follows:

- **path** = Path to parent directory of Fastq files
- **group** = Age group that sample belongs to e.g. _GROUP_
- **sample** = Sample identifier e.g. _TEST-0001_
- **filename** = Fastq file basename e.g. _TEST-0001-XT_
- **R1/R2** = Fastq filename (e.g `TEST-0001-XT_1.fastq.gz`)


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
```

