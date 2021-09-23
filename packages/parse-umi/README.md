# parse-umi

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
The command line script `parse-umi` is provided by this package for parsing UMI patterns from either read1 or read2 and labeling the read header with that information.


```
usage: parse-umi [-h] -1 READ1 -2 READ2 [-u1 UMI_READ1] [-u2 UMI_READ2]

Gather UMI information from read2

optional arguments:
  -h, --help            show this help message and exit
  -1 READ1, --read1 READ1
                        read1 fastq file
  -2 READ2, --read2 READ2
                        read2 fastq file
  -u1 UMI_READ1, --umi_read1 UMI_READ1
                        UMI Pattern for Read 1. for example XXNNNNNN (X =
                        nonUMI base, N=UMI bases)
  -u2 UMI_READ2, --umi_read2 UMI_READ2
                        UMI Pattern for Read 2. for example NNNNNNNNNNNN (X =
                        nonUMI base, N=UMI bases)						
```

### Example
Given the Example Below, Read labels will be modified to obtain the UMI sequences from both reads
```
parse-umi -1 ${read1} -2 ${read2} -u2 NNNNNNNNNNNNN -u1 XXXNNNNN

```
Read Label Example for the command above
```
@NB551443:36:HJ2FHAFXY:1:11101:5290:1029:UMI_R1__TTACC__32,32,36,36,36:UMI_R2__NTCTTCTATGTGT__2,32,32,32,32,36,36,36,36,36,36,36,36
```


#### Possible Read Labels
Read Labels have three potential formats, depending on the parameters passed

1) Both Read1 and Read2 UMI patterns `read_name:UMI1name__UMI1seq__UMI1qual:UMI2name__UMI2seq__UMI2qual`
1) Only Read1 UMI `read_name:UMI1name__UMI1seq__UMI1qual`
1) Only Read2 UMI `read_name:UMI1name__UMI1seq__UMI1qual:UMI2name__UMI2seq__UMI2qual`





