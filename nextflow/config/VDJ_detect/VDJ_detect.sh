#!/bin/bash
source activate Daedalus_env

VDJdetector -v $vSortBam -j $jSortBam -b ${sample} -i ${percentId} -r ${referenceData}


