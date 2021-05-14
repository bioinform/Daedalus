#!/bin/bash
source activate Daedalus_env

##run dedup with variable number of UMIs 
if [ "${umi_mode}" ==  "none" ]
then
    ipete-dedup -a $onTargetUMI -u ${umi_dist} -c ${cdr3_dist} -q ${min_qual} -l ${max_cdr3_len} -b ${bidding_ratio} -s ${max_steps} -p 1 -o $sample
fi

if [ "${umi_mode}" ==  "R1" ]
then
    ipete-dedup -a $onTargetUMI -u ${umi_dist} -c ${cdr3_dist} -q ${min_qual} -l ${max_cdr3_len} -b ${bidding_ratio} -s ${max_steps} -p 1 -o $sample -u1
fi
if [ "${umi_mode}" ==  "R2" ]
then
    ipete-dedup -a $onTargetUMI -u ${umi_dist} -c ${cdr3_dist} -q ${min_qual} -l ${max_cdr3_len} -b ${bidding_ratio} -s ${max_steps} -p 1 -o $sample -u2
fi
if [ "${umi_mode}" ==  "both" ]
then
    ipete-dedup -a $onTargetUMI -u ${umi_dist} -c ${cdr3_dist} -q ${min_qual} -l ${max_cdr3_len} -b ${bidding_ratio} -s ${max_steps} -p 1 -o $sample -u1 -u2
fi
    





