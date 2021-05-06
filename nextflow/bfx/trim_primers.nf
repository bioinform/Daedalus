nextflow.preview.dsl=2

include '../modules/trim_primers' params(params)

def trimPrimers(sample, primerRef, vPrimerAln, jPrimerAln) {
    (sample, trimmedFq, trimmedTsv, shortFq, shortTsv) =  trim_primers(sample, primerRef, vPrimerAln, jPrimerAln)
    return [sample, trimmedFq, trimmedTsv]    
}
