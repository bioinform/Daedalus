nextflow.preview.dsl=2

include '../modules/ipete_metrics' params(params)

def ipeteMetrics(sample, dedupReport) {
    (sample, metricsAll, metricsFunct) =  ipete_metrics(sample, dedupReport)
    return [sample, metricsAll, metricsFunct]    
}
