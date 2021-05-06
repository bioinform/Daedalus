nextflow.preview.dsl=2

finalDir = params.finalDir

process ipete_metrics {

    tag { sample }

    publishDir "${finalDir}/${sample}/ipete_metrics", mode: 'link'

    input:
	val(sample)
        file(dedupReport)
    
    output:
	val(sample)
        file("${sample}_cummulative_diversity_stats_all.csv")
        file("${sample}_cummulative_diversity_stats_functional.csv")
    
  script:
	template( task.ext.command )
      
}
