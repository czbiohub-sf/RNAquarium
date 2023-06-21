#!/bin/bash

working_dir="/hpc/scratch/group.theory/jparas/zf_pipeline"
sdir=${working_dir}/STAR_out

# reads lines from ${working_dir}/SRA_accession_list.1.27.23.txt, which is redirected after the done statement of this loop
while read -r accession; do
  # PE
  if [ -e ${sdir}/PE/${accession}/Log.final.out ] && [ -e ${sdir}/counts/${accession}/Log.final.out ] \
  && [ -e ${sdir}/PE/${accession}/Unmapped.out.mate1.gz ] && [ -e ${sdir}/PE/${accession}/Unmapped.out.mate2.gz ] \
  && [ -e ${sdir}/counts/${accession}/htseq-count.txt ]
  then
      htseq_size=$(wc -c ${sdir}/counts/${accession}/htseq-count.txt | awk '{print $1}')
      unmapped_1_size=$(gzip -c ${sdir}/PE/${accession}/Unmapped.out.mate1.gz | wc -c | awk '{print $1}')
      unmapped_2_size=$(gzip -c ${sdir}/PE/${accession}/Unmapped.out.mate2.gz | wc -c | awk '{print $1}')
      if  [ $htseq_size -gt 10000 ] && [ $unmapped_1_size -gt 100 ] && [ $unmapped_2_size -gt 100 ] && $(gzip -t ${sdir}/PE/${accession}/Unmapped.out.mate1.gz) && $(gzip -t ${sdir}/PE/${accession}/Unmapped.out.mate2.gz)
      then
          echo "${accession} STAR_out files all-OK" >> ${working_dir}/STARout.check.OK.txt # STAR_out OK, skip current job
      else
          echo "${accession} STAR_out files all-present but failed size/integrity checks. htseq-count.txt: ${htseq_size} Unmapped.out.mate1.gz: ${unmapped_1_size}" >> ${working_dir}/STARout.check.sizeFail.txt
      fi
  # SE
  elif [ -e ${sdir}/SE/${accession}/Log.final.out ] && [ -e ${sdir}/counts/${accession}/Log.final.out ] \
  && [ -e ${sdir}/counts/${accession}/htseq-count.txt ]  && [ -e ${sdir}/SE/${accession}/Unmapped.out.mate1.gz ]
  then
      htseq_size=$(wc -c ${sdir}/counts/${accession}/htseq-count.txt | awk '{print $1}')
      unmapped_1_size=$(gzip -c ${sdir}/SE/${accession}/Unmapped.out.mate1.gz | wc -c | awk '{print $1}')
      if [ $htseq_size -gt 10000 ] && [ $unmapped_1_size -gt 100 ] && $(gzip -t ${sdir}/SE/${accession}/Unmapped.out.mate1.gz)
      then
          echo "${accession} STAR_out files all-OK" >> ${working_dir}/STARout.check.OK.txt # STAR_out OK, skip current job
      else
          echo "${accession} STAR_out files all-present but failed size/integrity checks. htseq-count.txt: ${htseq_size} Unmapped.out.mate1.gz: ${unmapped_1_size}" >> ${working_dir}/STARout.check.sizeFail.txt
      fi
  fi
done < "${working_dir}/SRA_accession_list.1.27.23.txt"

