#!/bin/bash
source activate Daedalus_env
trim-primers -p $primerRef -v $vPrimerBam -j $jPrimerBam -b ${sample}

