# bam-streamer
A python package for reading pairs of alignments, combined from two separate bam files.



## Installation

## installation
```
pip install --upgrade git+http://ghe-rss.roche.com/pls-red-packages/bam-streamer.git
```


## Utilizing BamStreamer package

1) _Parse Bam File alignments_
read primary alignments, and high level stats, from a bam file

```python

from BamStreamer.BamUtils import read_alignments

for alignment in read_alignments("alignments.bam"):
    ref = alignment["reference_name"]
    start = alignment["start"]
    end = alignment["end"]
    #...

```

2) _Parse alignments from two bam files_

```python
from BamStreamer.BamPairs import read_two_bams

for alnPair = read_two_bams("aln1.bam", "aln1",  "aln2.bam", "aln2"):
    	ref1_name = alnPair["aln1"]["reference_name"]
	ref2_name = alnPair["aln2"]["reference_name"]
    	ref1_start = alnPair["aln1"]["start"]
	ref2_start = alnPair["aln2"]["end"]
				
```

