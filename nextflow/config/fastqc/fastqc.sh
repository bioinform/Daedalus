#!/bin/bash -e

# Load conda env within Docker
source /root/.bashrc || echo "Failed to source /root/.bashrc" >&2

    ${FASTQC_PATH}/fastqc ${read1} ${read2}
    unzip ${read1.name.replaceFirst(/.fastq$/, '_fastqc.zip')}
    unzip ${read2.name.replaceFirst(/.fastq$/, '_fastqc.zip')}

    mv ${read1.name.replaceFirst(/.fastq$/, '_fastqc')}/summary.txt \
        ${read1.name.replaceFirst(/.fastq$/, '_summary.txt')}
    mv ${read1.name.replaceFirst(/.fastq$/, '_fastqc')}/fastqc_data.txt \
        ${read1.name.replaceFirst(/.fastq$/, '_fastqc_data.txt')}

    mv ${read2.name.replaceFirst(/.fastq$/, '_fastqc')}/summary.txt \
        ${read2.name.replaceFirst(/.fastq$/, '_summary.txt')}
    mv ${read2.name.replaceFirst(/.fastq$/, '_fastqc')}/fastqc_data.txt \
        ${read2.name.replaceFirst(/.fastq$/, '_fastqc_data.txt')}
