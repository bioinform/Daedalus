#!/bin/bash -e

# Load conda env within Docker
source /root/.bashrc || echo "Failed to source /root/.bashrc" >&2

trim-primers -p $primerRef -v $vPrimerBam -j $jPrimerBam -b ${sample}
