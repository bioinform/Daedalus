# Daedalus

Nextflow pipeline for analysis of libraries prepared using the ImmunoPETE assay.

- [Daedalus](#daedalus)
  - [Install and Configure](#install-and-configure)
    - [Software Requirements](#software-requirements)
    - [Download git repo](#download-git-repo)
    - [Build Conda Environment (Optional)](#build-conda-environment-optional)
    - [Install](#install)
    - [Build Docker images](#build-docker-images)
    - [Configure images](#configure-images)
    - [Configure the pipeline](#configure-the-pipeline)
    - [Test Pipeline on a single sample](#test-pipeline-on-a-single-sample)
  - [Running Pipeline](#running-pipeline)
    - [Generate Manifest from Sample Sheet](#generate-manifest-from-sample-sheet)
    - [Submit Pipeline Run](#submit-pipeline-run)
    - [Output](#output)
  - [Workflow](#workflow)
  - [Methods](#methods)

## Install and Configure

Note... The Nextflow Config file must be configured for the queue.

### Software Requirements

- Python 3.6
- Java 8
- Nextflow 19.07.0, to run the pipeline
- UGE, for cluster job submission
- bats 0.4.0, for testing

### Download git repo

```bash
git clone git@github.com:bioinform/Daedalus.git
cd Daedalus
git checkout tags/${release-version}
```

### Build Conda Environment (Optional)

It's recommended to create a conda environment:

```bash
conda create -n Daedalus python=3.6
conda activate Daedalus
```

### Install

Within Daedalus directory, execute the following command.

```bash
pip install .
```

### Build Docker images

Due to license restriction, you will have to build the Bcl2fastq image using the Docker file.
Please refer to [Dockerhub](https://docs.docker.com/docker-hub/) for creating repo and pushing images.

```bash
docker build -t {dockerhub_username}/bcl2fastq:{version} -f Dockerfile_bcl2fastq .
docker push {dockerhub_username}/bcl2fastq:{version}
```

### Configure images

After building your own images, set the following params in the `nextflow/defaults-ipete.config` with your own images.

```javascript
params.bcl2fastq_docker = "{dockerhub_username}/bcl2fastq:{version}"
```

### Configure the pipeline

The pipeline runs on UGE cluster by default. If you install it on a different machine, modify the cluster settings in the `nextflow/nextflow.config` accordingly.

```javascript
ipete_docker {
    process.clusterOptions = { "-l h_vmem=${task.ext.vmem} -S /bin/bash -l docker_version=new -V" }
}
docker.runOptions = "-u=\$UID --rm -v /path/to/input_and_output:/path/to/input_and_output  -v /path/to/daedalus_repo:/path/to/daedalus_repo"
```

### Test Pipeline on a single sample
Once all the software has been installed and nextflow has been configured the pipeline bats test can be run. The bats test runs the pipeline on a single sample, from the paired fastq files provided:
- PBMC_1000ng_25ul_2_S6_R1_001.fastq.gz
- PBMC_1000ng_25ul_2_S6_R2_001.fastq.gz

In order the run the test, download both files from dropbox and move them into the data folder `Daedalus/data`. Once the data is available, run the test using the following commands:

```bash
cd test
bats single-sample-ipete.bats
```

An example of the pipeline output has also been provided: `PBMC_1000ng_25ul_2.tar.gz`

## Running Pipeline

Running the pipeline requires a complete flowcell worth of immunoPETE libraries.

### Generate Manifest from Sample Sheet

```bash
manifestGenerator = /path/to/Daedalus/pipeline_runner/manifest_generator.py
illuminaDir = /path/to/illumina/run_folder
sampleSheet = /path/to/sampleSheet.csv

python ${manifestGenerator} \
       --pipeline_run_id Daedalus_example_run \
       --sequencing_run_folder ${illuminaDir} \
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

python ${pipelineRunner} -g rssprbfprj --wait --resume -o ${outDir} Daedalus_example_manifest.csv
```

A `-g $group` needs to be provided to submit jobs to SGE cluster on SC1. 

### Output

At the specified output directory `${outDir}`, the analysis folder will be written using the `pipeline_run_id` "Daedalus_example_run"  

```bash
${outDir}/Daedalus_example_run
```

## Workflow

![workflow](docs/img/flowchart.png)

## Methods
Overview of the [Pipeline Methods](docs/Daedalus_methods.md) for key processing steps.
