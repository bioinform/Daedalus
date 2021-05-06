#!/usr/bin/env nextflow
nextflow.preview.dsl=2

//can also include workflow in separate segment
//include 'bfx/segment2-ipete' params(params)

include './bfx/subsample' params(params)
include './bfx/flash' params(params)
include './bfx/swifr_align' params(params)
include './bfx/parse_umi' params(params)
include './bfx/trim_primers' params(params)
include './bfx/VDJ_detect' params(params)
include './bfx/extract_umi' params(params)
include './bfx/ipete_dedup' params(params)
include './bfx/spikein_split' params(params)
include './bfx/ipete_reporter' params(params)

include ipete_metrics as ipete_metrics_native from './bfx/ipete_metrics' params(params)
include ipete_metrics as ipete_metrics_spikein from './bfx/ipete_metrics' params(params)

include swifr as swifrVprim from './modules/swifr_align' params(params)
include swifr as swifrJprim from './modules/swifr_align' params(params)
include swifr as swifrVgene from './modules/swifr_align' params(params)
include swifr as swifrJgene from './modules/swifr_align' params(params)
include swifr as swifrSpike from './modules/swifr_align' params(params)

//include one more module for pipeline switching

sample = params.sample
read1 = file(params.fastq1)
read2 = file(params.fastq2)

workflow {

    ///////////////////////////////////
    // quality filter, subsample reads
    ///////////////////////////////////
    //run seqtk
    (sample, subsampR1, subsampR2) = doSeqtk(sample, read1, read2)

    ////////////////////////////////////////////
    // label read headers with UMI (optionally)
    ////////////////////////////////////////////
    
    //If no UMI for read1 or read2, use "" instead of "NNNN..."
    if(params.umi_mode == true){
	(sample, read1UMI, read2UMI) = parseUMI(sample, subsampR1, subsampR2, params)
    } else{
	read1UMI = subsampR1
	read2UMI = subsampR2	
    }
    
    //combine reads with flash
    (sample, flashReads) = doFlash(sample, read1UMI, read2UMI)
    //if no flash? label read1 and read2 and process separately?    
    
    ///////////////
    //Trim Primers
    ///////////////
    if(params.trim_primers == true){
	// align V primers
	(sample2, vPrimerBam) = swifrVprim(sample, "align_primers", flashReads, file(params.vPrimerRef), "Vprimer", params.vPrimScore, params.vPrimKmer, params.vPrimFrac)

	// align J primers
	(sample, jPrimerBam) = swifrJprim(sample, "align_primers", flashReads, file(params.jPrimerRef), "Jprimer", params.jPrimScore, params.jPrimKmer, params.jPrimFrac)

	// trim primers
	(sample, trimmedFq, trimmedTsv) = trimPrimers(sample, file(params.primerRef), vPrimerBam, jPrimerBam)
    } else{
	trimmedFq = flashReads
    }
    
    /////////////////////////////////
    // align Spikeins (optionally)    
    /////////////////////////////////
    if(params.checkSpikein == true){
	(sample3, spikeBam) = swifrSpike(sample, "align_spikein", trimmedFq, file(params.spikeRef), "spikein", params.spikeScore, params.spikeKmer, params.spikeFrac)
    } 
    //////////////////////
    // align V & J genes
    //////////////////////

    // align V genes
    (sample2, vGeneBam) = swifrVgene(sample, "align_genes", trimmedFq, file(params.vGeneRef), "Vgene", params.vGeneScore, params.vGeneKmer, params.vGeneFrac)

    // align J genes
    (sample, jGeneBam) = swifrJgene(sample, "align_genes", trimmedFq, file(params.jGeneRef), "Jgene", params.jGeneScore, params.jGeneKmer, params.jGeneFrac)
    
    ///////////////////////////////
    // VDJ recombinant dectection
    ///////////////////////////////

    // VDJ detector
    (sample, onTargetTsv, onTargetFq, artifactTsv, artifactFq, alnSummary) = VDJdetect(sample, vGeneBam, jGeneBam, file(params.geneRef), params.minGeneIdentity)
                
    //////////////////////
    // dedup and consensus
    //////////////////////

    // extract UMI from read headers, add columns with UMI info
    // must be turned off without UMI
    if(params.umi_mode == true){
	(sample, onTargetUMI) = extractUMI(sample, onTargetTsv)	
    } else{
	onTargetUMI = onTargetTsv
    }

    // dedup types: 
    // 'CDR3' use CDR3 only
    // 'R1'- use Read1 UMI only
    // 'R2'- use Read2 UMI only
    // 'both'- use Read1 and Read2 UMIs
    if(params.umi_mode == true){
	(sample, dedupReport, onTargetDedup) = ipeteDedup(sample, onTargetUMI, params.umi_type)
    }
    
    ///////////////////////////////
    // split spikeins (optionally)    
    ///////////////////////////////
    if(params["checkSpikein"] == true){
	if(params.umi_mode == true){
	    (sample, dedupNative, dedupSpike, onTargetNative, onTargetSpike, vNativeBam, jNativeBam, vSpikeBam, jSpikeBam) = spikeinSplit(sample, spikeBam, vGeneBam, jGeneBam, dedupReport, onTargetDedup, file(params.geneRef))
	}	
    } else{	
	dedupNative = dedupReport
    }

    ///////////////////////////////////
    // CDR3 report & diversity metrics                            
    ///////////////////////////////////    
    if(params["checkSpikein"] == true){
	(sample, metricsAllSpike, metricsFunctSpike) =  ipete_metrics_spikein(sample, dedupSpike)
    }
    
    //(sample2, metricsAll, metricsFunct) =  ipete_metrics_native(sample, dedupNative)

    /////////////////////
    // filtered results
    /////////////////////    
    (sample, cdr3Report, diversityReport) = ipeteReporter(sample, dedupNative)    
    
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
