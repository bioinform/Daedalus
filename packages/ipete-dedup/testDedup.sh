
statFile=/sc1/groups/pls-redbfx/pipeline_runs/iPETE/production/2019-09-23-Expt56_AHJ2FHAFXY_Nextseq_immunoDb_replicate2/analysis/Exp56_PBMC_222_N714/filtered_reads/Exp56_PBMC_222_N714_on_target.tsv 

##statFile=/sc1/groups/pls-redbfx/immunoPETE/workflows/test/single-sample-ipete/testoutput/analysis/CCRF_CEM_10ng_N703/extract_umi/CCRF_CEM_10ng_N703_on_target_umi.tsv
sample=CCRF_CEM_10ng_N703

ipete-dedup -a ${statFile} -q 30 -o ${sample} --bidding_ratio 2 --steps 2 -u 1 -c 1 -l 100 -u2
