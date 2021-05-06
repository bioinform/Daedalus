#!/usr/bin/env nextflow


nextflow.preview.dsl=2

include '../modules/fastqc' params( params )

include ck_sample_metrics as ck_sample_metrics_fastqc1 from '../modules/ck_sample_metrics' params( params )
include gen_errorCodeSample as gen_errorCodeSample_fastqc1 from '../modules/ck_sample_metrics' params( params )
include kill_workflow as kill_workflow_fastqc1 from '../modules/ck_sample_metrics' params( params )

include ck_sample_metrics as ck_sample_metrics_fastqc2 from '../modules/ck_sample_metrics' params( params )
include gen_errorCodeSample as gen_errorCodeSample_fastqc2 from '../modules/ck_sample_metrics' params( params )
include kill_workflow as kill_workflow_fastqc2 from '../modules/ck_sample_metrics' params( params )
include fastqc as fastqc_1 from '../modules/fastqc' params( params )
include metricsFastqc as metricsFastqc_1 from '../modules/fastqc' params( params )


def doFastqcBeforeTrim(groupedDemux,preloadJson){
   metricsFastqc_1(fastqc_1( groupedDemux.merge(preloadJson).merge(Channel.from("before_trim")) )) 
}

def doFastqc(sampleInput, trimmed1, trimmed2, preloadJson){
   (fastqc_r1_qc_metrics, fastqc_r2_qc_metrics) = metricsFastqc( fastqc( sampleInput.merge(trimmed1).merge(trimmed2).merge(preloadJson).merge(Channel.from("after_trim"))) )

   val1 = ck_sample_metrics_fastqc1(sampleInput, fastqc_r1_qc_metrics, preloadJson )
   gen_errorCodeSample_fastqc1(sampleInput, val1)
   kill_workflow_fastqc1(val1)

   val2 = ck_sample_metrics_fastqc2(sampleInput, fastqc_r2_qc_metrics, preloadJson )
   gen_errorCodeSample_fastqc2(sampleInput, val2)
   kill_workflow_fastqc2(val2)
   
   result=[:]
   result['qcMetricsR1Fastqc']=fastqc_r1_qc_metrics
   result['qcMetricsR2Fastqc']=fastqc_r2_qc_metrics
   return result 
}

