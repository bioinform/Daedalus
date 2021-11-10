#!/bin/bash -e

# Load conda env within Docker
source /root/.bashrc || echo "Failed to source /root/.bashrc" >&2

##pipelineSummary=/sc1/groups/pls-redbfx/immunoPETE/develop/Daedalus/pipeline_runner/pipeline_summary.py

python ${params.pipeline_summary} -p ${analysis_dir}


