primer_file=/sc1/groups/pls-redbfx/users/dannebar/immunoPETE/analyses/20200203_Expt38_V2primers/V2.2_primers.csv
Valn=/sc1/groups/pls-redbfx/users/dannebar/immunoPETE/analyses/20200203_Expt38_V2primers/CDR3_analysis/Exp38_Strip1_1_Hut78_100ng_TRBVPool_10UMI_TRBJPool_0UMI.extendedFrags_vprimer_sort.bam
Jaln=/sc1/groups/pls-redbfx/users/dannebar/immunoPETE/analyses/20200203_Expt38_V2primers/CDR3_analysis/Exp38_Strip1_1_Hut78_100ng_TRBVPool_10UMI_TRBJPool_0UMI.extendedFrags_jprimer_sort.bam
sample=Exp38_Strip1_1_Hut78_100ng_TRBVPool_10UMI_TRBJPool_0UMI

python trimPrimers.py -p $primer_file -v $Valn -j $Jaln -b ${sample}
