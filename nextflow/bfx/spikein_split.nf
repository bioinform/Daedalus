nextflow.preview.dsl=2

include '../modules/spikein_split' params(params)

def spikeinSplit(sample, spikeBam, vSortBam, jSortBam, dedupReport, dedupReads, VJreference) {
    (sample, nativeDedup, spikeDedup, nativeReads, spikeReads, vNative, jNative, vSpike, jSpike) =  spikein_split(sample, spikeBam, vSortBam, jSortBam, dedupReport, dedupReads, VJreference)
    return [sample, nativeDedup, spikeDedup, nativeReads, spikeReads]    
}
