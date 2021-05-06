#!/bin/bash -e
bcl2fastq \
    --no-lane-splitting \
    -r ${task.cpus} \
    -p ${task.cpus} \
    -w ${task.cpus} \
    -R ${raw_dir} \
    -o . \
    --sample-sheet ${sample_sheet}
