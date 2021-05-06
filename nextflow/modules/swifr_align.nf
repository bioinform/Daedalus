nextflow.preview.dsl=2

finalDir = params.finalDir

process swifr {

    tag { sample }

    publishDir "${finalDir}/${sample}/${outdir}", mode: 'link'
    
    input:
	val(sample)
        val(outdir) 
        file(reads)        
        file(reference)
        val(outname)        
        val(alnScore)
        val(kmerSize)
        val(kmerFraction)
    
    output:
	val(sample)
        file("${outname}_sort.bam")

  script:
    template( task.ext.command )
}
