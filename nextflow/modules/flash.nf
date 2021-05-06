nextflow.preview.dsl=2

finalDir = params.finalDir

process flash {

    tag { sample }

    publishDir "${finalDir}/${sample}/flash", mode: 'link'

    input:
	val sample
        file read1
        file read2
    
    output:
	val(sample)
        file("${sample}.extendedFrags.fastq")
        file("${sample}.notCombined_1.fastq")
        file("${sample}.notCombined_2.fastq")
        file("${sample}.hist")

  script:
	template( task.ext.command )
      
}
