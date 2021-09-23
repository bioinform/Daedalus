# spikein-split 

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

## Usage
The command line script `spikein-split` is provided by this package for separating run results from spikein-molecules from
native molecules.

```
usage: spikein-split [-h] -s SPIKE_ALN -v V_ALN -j J_ALN -VJ VJ_UMI -d
                     DEDUP_REPORT -r REFERENCE_ANNOT -b BASENAME

gather amplicon stats from bam file

optional arguments:
  -h, --help            show this help message and exit
  -s SPIKE_ALN, --spike_aln SPIKE_ALN
                        Spike in alignment bam
  -v V_ALN, --v_aln V_ALN
                        V gene alignment bam
  -j J_ALN, --j_aln J_ALN
                        J gene alignment bam
  -VJ VJ_UMI, --VJ_umi VJ_UMI
                        V & J gene alignment summary file
  -d DEDUP_REPORT, --dedup_report DEDUP_REPORT
                        V & J gene alignment summary file
  -r REFERENCE_ANNOT, --reference_annot REFERENCE_ANNOT
                        reference annotation file
  -b BASENAME, --basename BASENAME
                        basename argument for report file name

```



