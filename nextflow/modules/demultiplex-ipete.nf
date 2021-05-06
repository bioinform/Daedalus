nextflow.preview.dsl=2

finalDir = params.finalDir

process bcl2fastq {
  tag { "bcl2fastq" }
  publishDir "${finalDir}/fastq", mode: 'link'

  input:
    file raw_dir
    file sample_sheet

  when:
    params.runBcl2Fastq

    output:
    file("*.fastq.gz")
    file("Reports")
    file("Stats")

  script:
    template( task.ext.command )
}
