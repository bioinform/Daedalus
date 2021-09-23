#dedupReport=/sc1/groups/pls-redbfx/users/dannebar/immunoPETE/analyses/20200203_Expt38_V2.2primers/nonTrim/Exp38_Strip1_5_PanT_100ng_TRBVPool_10UMI_TRBJPool_0UMI_cdr3_dedup_report.tsv

dedupReport=/sc1/groups/pls-redbfx/immunoPETE/workflows/test/single-sample-ipete/testoutput/analysis/Exp56_PBMC_222_N714/spikein_split/Exp56_PBMC_222_N714_cdr3_dedup_report_native.tsv

time ipete-reporter -f $dedupReport  -c 1 -s 2 -r 2 -b Exp56_PBMC_N714
