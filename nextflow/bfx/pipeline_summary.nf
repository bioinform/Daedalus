nextflow.preview.dsl=2

include '../modules/pipeline_summary' params(params)

def pipelineSummary(analysisDir) {
    sample_metrics = pipeline_summary(analysisDir)    
    return sample_metrics
}
