nextflow.preview.dsl=2

finalDir = params.finalDir

process pipeline_summary {

    tag { "sample_summary" }
    
    publishDir "${finalDir}/sample_summary", mode: 'link'

    input:
        file(analysis_dir)
    
    output:
        file("sample_metrics.csv")
    
  script:
	template( task.ext.command )
    
}
