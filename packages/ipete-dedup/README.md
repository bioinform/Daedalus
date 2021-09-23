# ipete-dedup

## Spack environment
This Package was built under spack version `v2.0.0`.
To load this spack version, edit your `.bashrc` to contain the following lines

```
APPS_DIR=/sc1/apps/spack/v2.0.0
APPS_MODULES="uge git miniconda3/4.3.14"
source /sc1/apps/setup/init.sh
```

After editing your `.bashrc` exit and re-login to update your spack version.

## Installation
1) Within your own conda environment, after cloning this repo
```
python setup.py install
```

2) Or, install provided conda enviroment
```
module load miniconda3
conda env create -f environment.yml
```

## Usage
The command line script `ipete-dedup` is provided by this package for deduplicating immunoPETE sequencing libraries. This parsing V and J gene alignments, finding successful rearrangements.

```
usage: ipete-dedup [-h] -a ALIGNMENT_SUMMARY -u UMI_DIST -c CDR3_DIST -q
                   MIN_QUAL -l MAX_CDR3_LENGTH -b BIDDING_RATIO -s STEPS
                   [-p PROCESSORS] -o OUT_BASENAME

Deduper designed for identifying UMI families and consensus sequences from
immunoPETE sequencing runs.

optional arguments:
  -h, --help            show this help message and exit
  -a ALIGNMENT_SUMMARY, --alignment_summary ALIGNMENT_SUMMARY
                        alignment summary file
  -u UMI_DIST, --umi_dist UMI_DIST
                        maximum edit distance allowed when searching UMIs
  -c CDR3_DIST, --cdr3_dist CDR3_DIST
                        maximum edit distance allowed when searching CDR3s
  -q MIN_QUAL, --min_qual MIN_QUAL
                        minimum average read quality to keep a read.
  -l MAX_CDR3_LENGTH, --max_cdr3_length MAX_CDR3_LENGTH
                        maximum CDR3 nt sequence length allowed (otherwise
                        read is filtered)
  -b BIDDING_RATIO, --bidding_ratio BIDDING_RATIO
                        The ratio when a smaller node will accept the bid of a
                        larger node.
  -s STEPS, --steps STEPS
                        number of edges to tranverse
  -p PROCESSORS, --processors PROCESSORS
                        number of processors used for deduping
  -o OUT_BASENAME, --out_basename OUT_BASENAME
                        basename for output summary reports: read report, cdr3
                        consensus report and umi-cdr3 node report
```

UMI and CDR3 sequences are used together to define UMI families for immunoPETE amplicons.
Three report files are written:
1) `read report`: per read UMI family information
2) `cdr3 cosensus report`: UMI family level information and consensus sequences
3) `umi-cdr3 node report`: a report of unique UMI-CDR3 combinations, used in the dedup process.


### Example
```
ipete-dedup -a ${on_target_reads_file} -q 30 -o ${sample_basename} --bidding_ratio 2 --steps 2 -u 1 -c 1 -l 100 

```



