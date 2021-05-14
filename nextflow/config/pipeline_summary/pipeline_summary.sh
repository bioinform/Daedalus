#!/bin/bash
source activate Daedalus_env

##pipeline-summary -p ${analysis_dir}

pipelineSummary=/sc1/groups/pls-redbfx/immunoPETE/develop/Daedalus/pipeline_runner/pipeline_summary.py

python \$pipelineSummary -p ${analysis_dir}


