#!/bin/bash
source activate Daedalus

VDJdetector -v $vSortBam -j $jSortBam -b ${sample} -i ${percentId} -r ${referenceData}


