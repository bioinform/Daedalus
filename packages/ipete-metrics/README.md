# ipete-metrics

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
Run the following within your own conda environment, after cloning this repo
```
python setup.py install
```

## Usage
The command line script `ipete-metrics` is provided by this package for calculating a variety of CDR3 diversity and entropy measures for UMI families.


```
usage: ipete-metrics [-h] -c CDR3_SUMMARY -b OUT_BASENAME

Gather cdr3 metrics from ipete stats file

optional arguments:
  -h, --help            show this help message and exit
  -c CDR3_SUMMARY, --cdr3_summary CDR3_SUMMARY
                        alignment summary file
  -b OUT_BASENAME, --out_basename OUT_BASENAME
                        output basename for report files

```

### Example
```
ipete-metrics -c ${cdr3Report} -b Exp56_PBMC_222_N714

```



