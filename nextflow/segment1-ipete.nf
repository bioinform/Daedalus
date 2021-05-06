#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include './bfx/demultiplex-ipete' params(params)

raw_dir = file(params.run_dir)
sample_sheet = file(params.sample_sheet)

workflow {
    
    fastqs = doBcl2Fastq(raw_dir, sample_sheet)
    
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
