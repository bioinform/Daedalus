# extract-umi 

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
The command line script `extract-umi` is provided by this package for parsing UMI labeled read headers. UMI string and quality values are extracted from read headers and added as new columns in the output tsv formatted table.

```
usage: extract-umi [-h] -r REPORTFILE -o OUTFILE

Parse UMI information from read header and store UMI information in new
columns

optional arguments:
  -h, --help            show this help message and exit
  -r REPORTFILE, --reportFile REPORTFILE
                        read report with UMI labeled headers
  -o OUTFILE, --outFile OUTFILE
                        output tsv file

```

### Example
```
extract-umi -r read_report.tsv -o read_report_umiExtract.tsv
```



