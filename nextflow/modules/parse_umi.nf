nextflow.preview.dsl=2

finalDir = params.finalDir

process parse_umi {

    tag { sample }

    publishDir "${finalDir}/${sample}/parse_umi", mode: 'link'

    input:
	val sample
        file read1
        file read2
        val umi_r1
        val umi_r2
    
    output:
	val(sample)
        file("${read1.name.replaceFirst(/.fastq$/,'')}_uid.fastq")
        file("${read2.name.replaceFirst(/.fastq$/,'')}_uid.fastq")

  script:
	template( task.ext.command )
      
}
