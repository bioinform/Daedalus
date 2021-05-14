#!/bin/bash
source activate Daedalus_env

spikein-split -s $spikeBam -v $vSortBam -j $jSortBam -VJ $dedupReads -d $dedupReport -r $VJreference -b $sample


