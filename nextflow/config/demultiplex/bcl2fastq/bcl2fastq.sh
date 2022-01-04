#!/bin/bash -e

# Load conda env within Docker
source /root/.bashrc || echo "Failed to source /root/.bashrc" >&2

bcl2fastq \
    --no-lane-splitting \
    -r ${task.cpus} \
    -p ${task.cpus} \
    -w ${task.cpus} \
    -R ${raw_dir} \
    -o . \
    --sample-sheet ${sample_sheet}
