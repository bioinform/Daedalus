nextflow.preview.dsl=2

include '../modules/swifr_align' params(params)

def doSwifr(sample, outdir, reads, reference, outname, alnScore, kmerSize, kmerFraction) {
    (sample, outFile) =  swifr(sample, outdir, reads, reference, outname, alnScore, kmerSize, kmerFraction)
    return [sample, outFile] 
}



