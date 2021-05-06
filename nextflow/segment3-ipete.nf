#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include './bfx/pipeline_summary' params(params)

analysisDir = file(params.analysis_dir)
workflow {
    
    sample_metrics = pipelineSummary(analysisDir)
    
} 


workflow.onComplete {
    println "Pipeline version: ${params.workflowVersion}"
    println "Pipeline started at: $workflow.start"
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Exit status: $workflow.exitStatus"
    println "Error message: $workflow.errorMessage"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    println "Exit status of the task that caused the workflow execution to fail: ${workflow.exitStatus}"
    println "Detailed error of the task that caused the workflow execution to fail: ${workflow.exitStatus}"
}
