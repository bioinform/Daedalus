# ipete-spikein-parser

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

2) Or, install the provided conda enviroment
```
module load miniconda3
conda env create -f environment.yml
```

## Usage
The command line script `ipeteSpikeins` is provided by this package with two primary functionalities.
1) `--action split`
By providing the argument `split` the scripts behavior is to separate spikein reads from sample reads. Both the V and J bam files will be split up into `native` and `spikein` bam files.

2) `--action count`
By providing the argument `count` a summary report of each spikein molecule counts is produced.


```
usage: ipeteSpikeins [-h] -s SPIKE_ALN [-v V_ALN] [-j J_ALN] [-VJ VJALN] -a
                     {split,count} -r REFERENCE_ANNOT -b BASENAME

Split spikein reads from native reads or count spikein molecules from a sample

optional arguments:
  -h, --help            show this help message and exit
  -s SPIKE_ALN, --spike_aln SPIKE_ALN
                        Spike in alignment bam
  -v V_ALN, --v_aln V_ALN
                        V gene alignment bam
  -j J_ALN, --j_aln J_ALN
                        J gene alignment bam
  -VJ VJALN, --VJaln VJALN
                        V & J gene alignment summary file
  -a {split,count}, --action {split,count}
                        fastq file
  -r REFERENCE_ANNOT, --reference_annot REFERENCE_ANNOT
                        reference annotation file
  -b BASENAME, --basename BASENAME
                        basename argument for report file name

```

### split example

```
ipeteSpikeins -s spike.bam -v v.bam -j j.bam -r immunoDB_reference.csv -b basename --action split

```

This `--action` produces several files using the same `basename`:
- `${basename}_Jgene_native.bam`
- `${basename}_Jgene_spikein.bam`
- `${basename}_Vgene_native.bam`
- `${basename}_Vgene_spikein.bam`


### count example

```
ipeteSpikeins -s $spikeAln --VJaln on_target_hits.csv -r immunoDB_reference.csv -b basename --action count

```

This `--action` produces a table of counts, using the `basename`:
- `${basename}_spikein_stats.csv`
