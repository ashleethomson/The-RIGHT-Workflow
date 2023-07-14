#!/usr/bin/env nextflow 

/*
################################################################################
Nextflow Definitions
################################################################################
*/

nextflow.enable.dsl=2
version = '0.0.1'

/*
################################################################################
Accessory functions to include
################################################################################
*/

// Adapt these functions (in lib/utilities.nf) to be part of the 'WorkflowMain' class
include {checkAndSetArgs; printArguments; getSampleFromCSV} from './lib/utilities.nf'

/*
################################################################################
Checking and printing input arguments
################################################################################
*/

WorkflowMain.callHelp(params.help, version) // Print help and remove default outdir ('./false')
checked_arg_map = checkAndSetArgs(params) // Check args for conflics
printArguments(checked_arg_map)           // Print pretty args to terminal

// Multiple assignment of read-regex and sample info to two variables
(reads, list) = getSampleFromCSV(checked_arg_map.samplesheet)

// Convert list of reads-regex to channel of reads
Channel
	.fromFilePairs(reads)
	.ifEmpty { exit 1, "Read files are empty"}
	.set { ch_reads }

// Convert list of sample information to channel
Channel
	.fromList(list)
	.set { ch_list }

// Channel of samples [ Filename, group, sample, path, [R1. R2] ]
ch_list.join(ch_reads).set { ch_samples }

/*
################################################################################
Implicit workflow: Run the QC-pipeline
################################################################################
*/

include { fusions } from './workflows/RIGHT-workflow.nf' params(checked_arg_map)

workflow {  
  fusions(
	  ch_samples
  )
}
