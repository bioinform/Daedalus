nextflow.preview.dsl=2

include '../modules/subsample' params(params)

def doSeqtk(sample, read1, read2) {
   return seqtk(sample, read1, read2)
}