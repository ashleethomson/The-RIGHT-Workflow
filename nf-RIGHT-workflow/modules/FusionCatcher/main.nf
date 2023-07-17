process fusioncatcher {

    tag { "fusioncatcher - ${filename}" } 
    publishDir "${outdir}/${group}/${filename}/fusioncatcher_v133", mode: 'copy'
    label 'process_fusioncatcher'

    input:
    tuple val(filename), val(group), val(sample), val(path), file(reads)
    val fusioncatcher_db
    val outdir
     
    output:
    file "${filename}_fusioncatcher.txt"
    file "${filename}_fusioncatcher_hg38.txt"
    file "${filename}_fusioncatcher_sequences.txt.zip"

    script:
    """
	fusioncatcher \
        -p ${task.cpus} \
        -d ${fusioncatcher_db} \
        -i ${reads[0]},${reads[1]} \
        -o \${PWD} \
        --skip-blat
			
	mv final-list_candidate-fusion-genes.hg19.txt ${filename}_fusioncatcher.txt
	mv final-list_candidate-fusion-genes.txt ${filename}_fusioncatcher_hg38.txt
	mv final-list_candidate-fusion-genes_sequences.txt.zip ${filename}_fusioncatcher_sequences.txt.zip
    """
}
