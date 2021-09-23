# trim-primers

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
Within your own conda environment, after cloning this repo
```
python setup.py install
```

## Usage
The command line script `trim-primers` is provided by this package for taking pairs of primer alignments and trimming reads of those alignments.


```
usage: trim-primers [-h] -p PRIMER_INFO -b BASENAME -v V_ALN -j J_ALN

gather amplicon stats from bam file

optional arguments:
  -h, --help            show this help message and exit
  -p PRIMER_INFO, --primer_info PRIMER_INFO
                        primer information file
  -b BASENAME, --basename BASENAME
                        basename for output files
  -v V_ALN, --v_aln V_ALN
                        V probe alignment bam file
  -j J_ALN, --j_aln J_ALN
                        J probe alignment bam file

```

### Example
```
trim-primers -p $primer_reference_csv -v $Valn -j $Jaln -b $basename

```



