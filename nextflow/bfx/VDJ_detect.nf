nextflow.preview.dsl=2

include '../modules/VDJ_detect' params(params)

def VDJdetect(sample, vSortBam, jSortBam, referenceData, percentId) {
    (sample, onTargetTsv, onTargetFq, artifactTsv, artifactFq, alnSummary) =  VDJ_detect(sample, vSortBam, jSortBam, referenceData, percentId)
    return [sample, onTargetTsv, onTargetFq, artifactTsv, artifactFq, alnSummary]    
}
