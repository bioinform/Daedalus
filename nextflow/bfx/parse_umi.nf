nextflow.preview.dsl=2

include '../modules/parse_umi' params(params)

def parseUMI(sample, read1, read2, params) {
    (sample, read1UMI, read2UMI) =  parse_umi(sample, read1, read2, params.umi1, params.umi2)
    return [sample, read1UMI, read2UMI]    
}
