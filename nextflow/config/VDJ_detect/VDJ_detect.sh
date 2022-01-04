#!/bin/bash -e

# Load conda env within Docker
source /root/.bashrc || echo "Failed to source /root/.bashrc" >&2

VDJdetector -v $vSortBam -j $jSortBam -b ${sample} -i ${percentId} -r ${referenceData}


