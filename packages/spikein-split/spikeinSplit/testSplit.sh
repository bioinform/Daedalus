##/sc1/groups/pls-redbfx/immunoPETE/workflows/test/single-sample-ipete/testoutput/analysis/CCRF_CEM_10ng_N703

dir=..
ref=/sc1/groups/pls-redbfx/pipeline_runs/iPETE/production/2019-09-23-Expt60_shiny_test/reference/immunoDB_VJtargets.csv
spike_aln=${dir}/spikein_sort.bam
v_aln=${dir}/Vgene_sort.bam
j_aln=${dir}/Jgene_sort.bam
UMIReport=${dir}/Exp56_PBMC_222_N714_on_target_umi_groups.tsv
DedupReport=${dir}/Exp56_PBMC_222_N714_cdr3_dedup_report.tsv
python spikeinSplit.py -s $spike_aln -v $v_aln -j $j_aln -VJ $UMIReport -d ${DedupReport} -r $ref -b test_split_out
