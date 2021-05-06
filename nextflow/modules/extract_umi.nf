nextflow.preview.dsl=2

finalDir = params.finalDir

process extract_umi {

    tag { sample }

    publishDir "${finalDir}/${sample}/extract_umi", mode: 'link'

    input:
	val(sample)
        file(readReport)
    
    output:
	val(sample)
        file("${sample}_on_target_umi.tsv")
    
  script:
	template( task.ext.command )
      
}
