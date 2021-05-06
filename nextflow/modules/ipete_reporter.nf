nextflow.preview.dsl=2

finalDir = params.finalDir

//params["cdr3_edit_dist"]
//params["bidding_ratio"]
//params["max_steps"]

process ipete_reporter {

    tag { sample }

    publishDir "${finalDir}/${sample}/ipete_reporter", mode: 'link'

    input:
	val(sample)
        file(dedupReport)
        val(cdr3_edit_dist)
        val(max_steps)
        val(bidding_ratio)
    
    output:
	val(sample)
        file("${sample}_functional_D99umi_CDR3_report.tsv")
        file("${sample}_CDR3_diversity_stats.tsv")
    
  script:
	template( task.ext.command )
      
}
