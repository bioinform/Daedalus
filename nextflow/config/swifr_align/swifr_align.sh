#!/bin/bash -e

${params.swifr} -v -f $reads -q $reference -k $kmerSize -m 5 -s $alnScore -p ${task.cpus} -F $kmerFraction -o ${outname}
samtools view -Sb ${outname}.sam > ${outname}.bam
samtools sort -n -@ ${task.cpus} ${outname}.bam ${outname}_sort
sleep 10

#echo `date +"%T"` > mytime.txt
#sleep 20 
#echo `date +"%T"` > mytimeEnd.txt
#touch ${outname}_sort.bam
