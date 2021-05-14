#!/bin/bash
source activate Daedalus_env

if [ -z "${umi_r1}" ]
then
      parse-umi -1 ${read1} -2 ${read2} -u2 ${umi_r2}
else
    ##if no umi on R2, check R1 only
    if [ -z "${umi_r2}" ]
    then
	parse-umi -1 ${read1} -2 ${read2} -u1 ${umi_r1} 	
    else
	##otherwise, check both R1 and R2 for UMI
	parse-umi -1 ${read1} -2 ${read2} -u1 ${umi_r1} -u2 ${umi_r2}
    fi    
fi

