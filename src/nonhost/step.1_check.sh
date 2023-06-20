#!/bin/bash

working_dir="/hpc/scratch/group.theory/jparas/zf_pipeline"
sdir=${working_dir}/STAR_out

# PE
if [ -e ${sdir}/PE/${accessions[$idx]}/Log.final.out ] && [ -e ${sdir}/counts/${accessions[$idx]}/Log.final.out ] \
&& [ -e ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate1.gz ] && [ -e ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate2.gz ] \
&& [ -e ${sdir}/counts/${accessions[$idx]}/htseq-count.txt ]
then 
    countFileSize=$(wc -c ${sdir}/counts/${accessions[$idx]}/htseq-count.txt | awk '{print $1}')
    countFileSize2=$(gzip -c ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate1.gz | wc -c | awk '{print $1}')
    if  [ $countFileSize -gt 100 ] && [ $countFileSize2 -gt 100 ] && $(gzip -t ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate1.gz) && $(gzip -t ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate2.gz)
    then
        echo "${accessions[$idx]} STAR_out files all-OK" >> ${working_dir}/STARout.check.OK.txt # STAR_out OK, skip current job
        exit 0
    else
        echo "${accessions[$idx]} STAR_out files all-present but failed size/integrity checks. htseq-count.txt: ${countFileSize} Unmapped.out.mate1.gz: ${countFileSize2}" >> ${working_dir}/STARout.check.sizeFail.txt
        exit 0
    fi
fi
# SE
if [ -e ${sdir}/SE/${accessions[$idx]}/Log.final.out ] && [ -e ${sdir}/counts/${accessions[$idx]}/Log.final.out ] \
&& [ -e ${sdir}/counts/${accessions[$idx]}/htseq-count.txt ]  && [ -e ${sdir}/SE/${accessions[$idx]}/Unmapped.out.mate1.gz ] 
then  
    countFileSize=$(wc -c ${sdir}/counts/${accessions[$idx]}/htseq-count.txt | awk '{print $1}')
    countFileSize2=$(gzip -c ${sdir}/SE/${accessions[$idx]}/Unmapped.out.mate1.gz | wc -c | awk '{print $1}')
    if [ $countFileSize -gt 100 ] && [ $countFileSize2 -gt 100 ] && $(gzip -t ${sdir}/SE/${accessions[$idx]}/Unmapped.out.mate1.gz)
    then
        echo "${accessions[$idx]} STAR_out files all-OK" >> ${working_dir}/STARout.check.OK.txt # STAR_out OK, skip current job
        exit 0
    else
        echo "${accessions[$idx]} STAR_out files all-present but failed size/integrity checks. htseq-count.txt: ${countFileSize} Unmapped.out.mate1.gz: ${countFileSize2}" >> ${working_dir}/STARout.check.sizeFail.txt
        exit 0
    fi
fi