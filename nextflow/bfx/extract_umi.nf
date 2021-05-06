nextflow.preview.dsl=2

include '../modules/extract_umi' params(params)

def extractUMI(sample, readReport) {
    (sample, outFile) =  extract_umi(sample, readReport)
    return [sample, outFile]
}
