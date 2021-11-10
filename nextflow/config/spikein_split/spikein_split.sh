#!/bin/bash -e

# Load conda env within Docker
source /root/.bashrc || echo "Failed to source /root/.bashrc" >&2

spikein-split -s $spikeBam -v $vSortBam -j $jSortBam -VJ $dedupReads -d $dedupReport -r $VJreference -b $sample


