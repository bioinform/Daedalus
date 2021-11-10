#!/bin/bash -e

# Load conda env within Docker
source /root/.bashrc || echo "Failed to source /root/.bashrc" >&2

ipete-reporter -f $dedupReport -b $sample -c $cdr3_edit_dist -s $max_steps -r $bidding_ratio



