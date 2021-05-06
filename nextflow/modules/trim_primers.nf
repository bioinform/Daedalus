nextflow.preview.dsl=2

finalDir = params.finalDir

process trim_primers {

    tag { sample }

    publishDir "${finalDir}/${sample}/trim_primers", mode: 'link'

    input:
	val sample
        file primerRef
        file vPrimerBam
        file jPrimerBam
    
    output:
	val(sample)
        file("${sample}_trim_primers.fastq")
        file("${sample}_trim_primers.tsv")
        file("${sample}_too_short.fastq")
        file("${sample}_too_short.tsv")

  script:
	template( task.ext.command )
      
}
