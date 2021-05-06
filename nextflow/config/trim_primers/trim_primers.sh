#!/bin/bash
source activate Daedalus
trim-primers -p $primerRef -v $vPrimerBam -j $jPrimerBam -b ${sample}

