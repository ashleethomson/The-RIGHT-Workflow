/*
RNA-seq fusion calling with new programs
	- FASTQ to FusionCatcher, Arriba and STAR-Fusion
*/

// Import processes used by sub-workflow
include { arriba } from './../modules/arriba/2.1.0/Arriba'
include { starFusion } from './../modules/STARFusion/1.10.0/star-fusion'
include { fusioncatcher_v133 } from './../modules/fusioncatcher/1.33/fusioncatcher'

/*
Call the sub-workflow: fusions
*/

workflow fusions {
	take:
		ch_samples
	
	main:

		starFusion(
			ch_samples,
			params.starfuse_ctat,
			params.outdir
		)

		arriba(
			ch_samples,
			params.assembly,
			params.gtf,
			params.arriba_blacklist,
			params.arriba_knownfus,
			params.arriba_proteindom,
			params.arriba_staridx,
			params.outdir
		)

		fusioncatcher_v133(
			ch_samples
			params.fusioncatcher_db
			params.outdir
		)
}
