nextflow.preview.dsl=2

finalDir = params.finalDir

process seqtk {
  tag { sample }
  publishDir "${finalDir}/${sample}/seq/seqtk", mode: 'link'

  input:
    val sample
    file read1
    file read2

  output:
    val sample
    file ( "${read1.name.replaceFirst(/.fastq.*$/, '')}_seqtk.fastq" )
    file ( "${read2.name.replaceFirst(/.fastq.*$/, '')}_seqtk.fastq" )

  script:
    template( task.ext.command )
}