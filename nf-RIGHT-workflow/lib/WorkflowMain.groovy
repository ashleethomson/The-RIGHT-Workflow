class WorkflowMain {
    public static String printVersion(String version) {
		def message = ''
        message += """
		==============================================================
			   ALL GENE-FUSION NEXTFLOW PIPELINE ${version}		
		==============================================================
		""".stripIndent()

        println(message)
    }

	public static String printHelpMessage() {
		def message = ''
		message += """
		Nextflow Arguments:
			-profile <str>			Which Nextflow profile to use: Should ALWAYS be 'conda,slurm'
			-N <str>			    Email that a notification of completion will be sent to
		
		Mandatory Arguments:
			--outdir <str>			Path to output directory
			--samplesheet <str>		Path to sample sheet
			--email <str>			Your email
			--genome <str>          Path to GRCh37 genome
			--annotaion <str>       Path to GRCh37 annotaion
			--workDir <str>         Path to work directory location
		""".stripIndent()
        println(message)
	}

    public static void callHelp(boolean help, String version) {
		if (help) {
            def temp_outdir = new File('false')
            if (temp_outdir.exists()) {
                temp_outdir.deleteDir()
            }
			printVersion(version)
			printHelpMessage()
			System.exit(0)
		}
	}
}