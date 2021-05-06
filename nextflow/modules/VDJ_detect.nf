nextflow.preview.dsl=2

finalDir = params.finalDir

process VDJ_detect {

    tag { sample }

    publishDir "${finalDir}/${sample}/VDJ_detect", mode: 'link'

    input:
	val sample
        file vSortBam
        file jSortBam
        file referenceData
        val percentId
    
    output:
	val(sample)
        file("${sample}_on_target.tsv")
        file("${sample}_on_target.fastq")
        file("${sample}_artifacts.tsv")
        file("${sample}_artifacts.fastq")
        file("${sample}_alignment_summary.csv")
        
  script:
	template( task.ext.command )
      
}
