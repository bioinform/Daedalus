nextflow.preview.dsl=2

include '../modules/flash' params(params)

def doFlash(sample, read1, read2) {
    (sample, flashReads, noFlash1, noFlash2, flashHist) =  flash(sample, read1, read2)
    return [sample, flashReads]    
}
