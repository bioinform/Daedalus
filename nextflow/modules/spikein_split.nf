nextflow.preview.dsl=2

finalDir = params.finalDir

process spikein_split {

    tag { sample }

    publishDir "${finalDir}/${sample}/spikein_split", mode: 'link'

    input:
	val(sample)
        file(spikeBam)
        file(vSortBam)
        file(jSortBam)
        file(dedupReport)
        file(dedupReads)
        file(VJreference)
    
    output:
	val(sample)
        file("${sample}_cdr3_dedup_report_native.tsv")
        file("${sample}_cdr3_dedup_report_spikein.tsv")
        file("${sample}_on_target_umi_groups_native.tsv")
        file("${sample}_on_target_umi_groups_spikein.tsv")
        file("${sample}_Vgene_native.bam")
        file("${sample}_Jgene_native.bam")
        file("${sample}_Vgene_spikein.bam")
        file("${sample}_Jgene_spikein.bam")    
    
  script:
	template( task.ext.command )
      
}
