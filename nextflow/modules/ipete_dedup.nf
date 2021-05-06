nextflow.preview.dsl=2

finalDir = params.finalDir

process ipete_dedup {

    tag { sample }

    publishDir "${finalDir}/${sample}/dedup", mode: 'link'

    input:
	val sample
        file onTargetUMI
        val umi_mode
        val max_cdr3_len
        val max_steps
        val bidding_ratio
        val umi_dist
        val cdr3_dist
        val min_qual
    output:
	val(sample)
        file("${sample}_cdr3_dedup_report.tsv")
        file("${sample}_on_target_umi_groups.tsv")
    
  script:
	template( task.ext.command )
      
}
