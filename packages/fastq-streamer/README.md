# fastq-streamer

## installation
```
pip install --upgrade git+http://ghe-rss.roche.com/pls-red-packages/fastq-streamer.git
```

## Utilizing the FastqReader class

1) _Import FastqReader class_

```python
from FastqStreamer.FastqReader import FastqReader
```


2) _Parse Fastq Files_
_read each line of the fastq file as a dictionary of read information)

```python

FqReader = FastqReader()
filter_out = FqReader.create_fastq_handle("filtered_reads.fastq")

for read in FqReader.iter_reads("reads.fastq"):
	read_name = read["read_name"]
	seq = read["seq"]
	qual =  read["qual"]
	## filter reads: quality, length, barcode, etc...
	FqReader.write_fastq(read, filter_out)

##close fastq file handles
FqReader.close_file_handles()

```


3) _Parse paired end Fastq files_

```python
FqReader = FastqReader()

filter1_out = FqReader.create_fastq_handle("filtered_read1.fastq")
filter2_out = FqReader.create_fastq_handle("filtered_read2.fastq")

for read in FqReader.iter_pairs("read1.fastq", "read2.fastq")
    	read_name = read["read1"]["read_name"]
		seq_r1 = read["read1"]["seq"]
		qual_r1 =  read["read1"]["qual"]
		seq_r2 = read["read2"]["seq"]
		qual_r2 =  read["read2"]["qual"]
		## filter reads ...
		FqReader.write_fastq(read["read1"], filter1_out)
		FqReader.write_fastq(read["read2"], filter2_out)
		
FqReader.close_file_handles()
		
```



