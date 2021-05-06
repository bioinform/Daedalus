## Pipeline options

- checkSpikeIn: `true or false`
Aligns spike-in barcode sequences against reads, all spike-in reads detected
are analyzed separately in the pipeline (default = true)

- subsample: `int`
By default, all reads are analyzed by the pipeline (subsample = 1). Setting this parameter will subset
the sample to the desired number of read pairs, useful for debugging and quickly checking sample quality.

- numProcessors: `int`
number of processors to use for parallelized steps.

- minAlnScore: `int`
The minimum alignment score threshold to report an alignment. Recommended > 16

- uid_length: `int`
The uid sequence length on the vgene primer. (default = 13)

- trimreads:`true or false`
Trim reads of adapter sequence and by read quality using trimmomatic. 

- trim_qual: `int`
The phred quality-score threshold used by trimmomatic to trim reads.

- runFastqc: `true or false`
Run fastqc on the samples (default = true)

- runPresto: `true or false`
Find consensus using PRESTO. (default = true)

- runCobbSalad: `true or false`
Find consensus using CobbSalad. (default = false)

