# VDJ-recombinant-detection

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
The command line script `VDJdetector` is provided by this package for parsing V and J gene alignments, finding successful rearrangements.

```
usage: VDJdetector [-h] -r REFERENCE_INFO -b BASENAME -v VGENE_ALN -j
                   JGENE_ALN [-i IDENTITY]

Identify VDJ recombinants from V and J gene alignments

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE_INFO, --reference_info REFERENCE_INFO
                        immunoDB reference annotation
  -b BASENAME, --basename BASENAME
                        basename for output files
  -v VGENE_ALN, --vgene_aln VGENE_ALN
                        V gene alignment bam file
  -j JGENE_ALN, --jgene_aln JGENE_ALN
                        J gene alignment bam file
  -i IDENTITY, --identity IDENTITY
                        minimum gene identity to be considered a good
                        alignment

```

CDR3 sequences can be discovered from V and J gene segment alignments that are performed on the same set of reads. V and J segments that are aligned must also have CDR3 position annotation, which are available in the immunoDB reference (http://ghe-rss.roche.com/plsRED-Bioinformatics/immunoDB). `VDJdetector` utilizes the V and J gene bam files along with the reference annotation to find VDJ recombinants, CDR3 sequence and quality, and D-regions. `VDJdetector` also requires both `vgene` and `jgene` bam files to be sorted by read name (`samtools sort -n`) in order to enable alignments to the same read to be parsed together.

### Example
```
VDJdetector -v vgenes.bam -j jgenes.bam -b sample_basename --identity 95 -r immunoDB_reference.csv
```



