# Daedalus

Nextflow pipeline for analysis of libraries prepared using the ImmunoPETE assay.

- [Daedalus](#daedalus)
    - [Download](#download)
    - [Build Conda Environment](#build-conda-environment)
    - [Test Pipeline](#test-pipeline)
    - [Run Pipeline](#run-pipeline)
        - [Load Environment](#load-the-environment)
        - [Generate Manifest from Sample Sheet](#generate-manifest-from-sample-sheet)
        - [Submit Pipeline Run](#submit-pipeline-run)
        - [Output](#output)
    - [Workflow](#workflow)
    - [Methods](#methods)

Note... The Nextflow Config file must be configured for the queue.

## Software Requirements 
- built on a linux server: CentOS Linux release 7.7.1908 (Core)
- miniconda3, for package management
- nextflow 19.07.0, to run the pipeline 
- uge, for cluster job submission

## Install SWIFR aligner
A smith waterman alignment implemention (c++) was developed and is used to identfy primers and V/J gene segements
Full README for ![swifr](packages/swifr/README.md) and how to install.

## Download git rep
```bash
git clone git@github.com:bioinform/Daedalus.git
cd Daedalus
```

## Install SWIFR aligner
Full README for ![swifr](packages/swifr/README.md) and how to install.


## Build Conda Environment

Build the conda environment for running the pipeline:
```bash
conda env create -f environment.yml
```

##install python packages in the loaded conda ENV
must be performed within the cloned directory

```bash
conda activate Daedalus
./install_packages.sh
```

## Nextflow configuration
Nextflow must be configured for each systm. The ipete profile in the nextflow config file `./nextflow/nextflow.config` should be updated accordingly.


## Test Pipeline on a single sample

```bash
module load bats
module load nextflow_latest/19.07.0
cd test
bats single-sample-ipete.bats
```

## Running Pipeline

### Load the Environment

Activate the conda environment and be sure nextflow is loaded  

```bash
conda activate Daedalus
```

### Generate Manifest for ImmunoPETE Run from the Sample Sheet

```bash
manifestGenerator = /path/to/Daedalus/pipeline_runner/manifest_generator.py
illuminaDir = /path/to/illumina/run_folder
sampleSheet = /path/to/sampleSheet.csv

python ${manifestGenerator} \
       --pipeline_run_id Daedalus_example_run \
       --sequencing_run_folder ${illuminaDir} \
       --sequencing_platform NextSeq \
       --output Daedalus_example_manifest.csv \
       --subsample 1 \
       --umi_mode True \
       --umi2 'NNNNNNNNN' \
       --umi_type R2 \
	   ${sampleSheet}

```

The manifest file contains all parameters needed for the pipeline to run. Sample specific tuning of parameters or any updates to the parameters can be acheived by editing the manifest file generated. After edits are complete, the pipeline can be submitted using the manifest file alone.

### Submit Pipeline Run

Using the output from Manifest Generator `Daedalus_example_manifest.csv` pipeline runs can be submitted using the script: pipeline_runner.py.

```bash
pipelineRunner=/path/to/Daedalus/pipeline_runner/pipeline_runner.py
outDir=/path/to/analysis/output

python ${pipelineRunner} -g rssprbf --wait --resume -o ${outDir} Daedalus_example_manifest.csv
``````

A `-g $group` needs to be provided to submit jobs to SGE cluster on SC1. 

### Output

At the specified output directory `${outDir}`, the analysis folder will be written using the `pipeline_run_id` "Daedalus_example_run"  

```bash
${outDir}/Daedalus_example_run
```

## Workflow

![workflow](docs/img/flowchart.png)








