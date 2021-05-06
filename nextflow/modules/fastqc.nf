nextflow.preview.dsl=2

ctx = params.ctx
finalDir = ctx.getFinalDir() 
params['fastqcVersion'] = ctx.getFastqcVersion()
globals= ctx.getGlobals()
sampleData= ctx.getSampleData()
params['oncoNgsUtilsVersion'] = ctx.getOncoNgsUtilsVersion()
params['rebelJVMArgs'] = ctx.getRebelJVMArgs()

process fastqc {
  tag { "${sample}-${order}" }
  publishDir "${finalDir}/${sample}/fastqc/${order}", mode: 'link'

  afterScript 'sync; sync; sync'

  input:
    set val(sample), file(read1), file(read2), file(preload), val(order)

  output:
    val(sample)
    file ('*.html')
    file ('*.zip')
    file ( "${read1.name.replaceFirst(/.fastq/, '')}_summary.txt" ) 
    file ( "${read1.name.replaceFirst(/.fastq/, '')}_fastqc_data.txt" )
    file ( "${read2.name.replaceFirst(/.fastq/, '')}_summary.txt" )
    file ( "${read2.name.replaceFirst(/.fastq/, '')}_fastqc_data.txt" )
    val(order)

  script:
      template( task.ext.command )
}

process metricsFastqc {
  tag { "${sample}_${order}" }
  publishDir "${finalDir}/${sample}/fastqc/${order}", mode: 'link'

  input:
    val(sample)
    file(html)
    file(zip)
    file(r1_fastqc_summary_file)
    file(r1_fastqc_data_file)
    file(r2_fastqc_summary_file)
    file(r2_fastqc_data_file)
    val(order)

  output:
    file("${sample}_fastqc_r1_qc_metrics.txt") //into qcMetricsR1Fastqc
    file("${sample}_fastqc_r2_qc_metrics.txt") // into qcMetricsR2Fastqc

  script:
      template( task.ext.command )
}

