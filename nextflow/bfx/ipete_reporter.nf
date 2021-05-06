nextflow.preview.dsl=2

include '../modules/ipete_reporter' params(params)

def ipeteReporter(sample, dedupReport) {
    (sample, cdr3Report, diversityReport) =  ipete_reporter(sample, dedupReport, params.cdr3_edit_dist, params.max_steps, params.bidding_ratio)
    return [sample, cdr3Report, diversityReport]    
}
