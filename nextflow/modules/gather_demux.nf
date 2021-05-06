nextflow.preview.dsl=2

finalDir = params.finalDir

process gather_demux {
  tag { "bcl2fastq" }
  publishDir "${finalDir}/fastq", mode: 'link'

  input:
    file raw_dir
    file sample_sheet

  when:
    params.runBcl2Fastq

  output:
    file "**/*.fastq.gz"
    file "Undetermined*.fastq.gz"
    file "Reports"
    file "Stats"

  script:
    template( task.ext.command )
}
