#!/bin/bash
flash ${read1} ${read2} -O -o ${sample} -M 250
##and uncombined read1
cat ${sample}.notCombined_1.fastq >> ${sample}.extendedFrags.fastq



