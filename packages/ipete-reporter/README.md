# ipete-reporter

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
The command line script `ipete-reporter` is provided by this package for filtering the immunoPETE pipeline output, producing high quality UMI families, CDR3 clusters  and CDR3 diversity statistics for functional T/B-cell clones

```
usage: ipete-reporter [-h] -f CDR3_REPORT -c CDR3_DIST -s MAX_STEPS -r
                      BIDDING_RATIO -b BASENAME

Filter CDR3 dedup report (D99 + functional), reporting high quality CDR3, CDR3
clusters, and divserity stats

optional arguments:
  -h, --help            show this help message and exit
  -f CDR3_REPORT, --cdr3_report CDR3_REPORT
                        cdr3 dedup report
  -c CDR3_DIST, --cdr3_dist CDR3_DIST
                        cdr3 clustering, edit distance
  -s MAX_STEPS, --max_steps MAX_STEPS
                        cdr3 clustering, scoring steps/distance
  -r BIDDING_RATIO, --bidding_ratio BIDDING_RATIO
                        cdr3 clustering, bidding ratio
  -b BASENAME, --basename BASENAME
                        output basename for report files

```

### Example
```
ipete-reporter -f ${cdr3Report} -b Exp56_PBMC_222_N714 -c 1 -s 2 -r 2

```



