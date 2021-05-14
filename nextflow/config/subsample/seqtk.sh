#!/bin/bash
# Default params.subsample should be set to 1 so seqtk won't do subsampling.
source activate Daedalus_env

if [ ${params.subsample} == 1 ] 
then
    gunzip -c $read1 > read1.fastq
    gunzip -c $read2 > read2.fastq 
    mv read1.fastq ${read1.name.replaceFirst(/.fastq.*$/, '_seqtk.fastq')}
    mv read2.fastq ${read2.name.replaceFirst(/.fastq.*$/, '_seqtk.fastq')}
else
    seqtk sample -s ${params.subsampleSeed} ${read1} ${params.subsample} > \
	  ${read1.name.replaceFirst(/.fastq.*$/, '_seqtk.fastq')}
    seqtk sample -s ${params.subsampleSeed} ${read2} ${params.subsample} > \
	  ${read2.name.replaceFirst(/.fastq.*$/, '_seqtk.fastq')}

fi
