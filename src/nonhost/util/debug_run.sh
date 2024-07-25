#!/bin/bash
sed -nr 's/.*\(([^)]+).*terminated with an error exit status.*/\1/p' $1 | sort | uniq |
while read -r accession; do
# read process name, accession, slurm job id, working directory for given accession
echo "ACCESSION: $accession" | tee -a "$accession.log"
grep -o "[^ ]\+ ($accession) > jobId: [0-9]\+; workDir: [^ ]\+" $1 | cut -d' ' -f1,2,5,7 | tr -d ";" |
  while read -r proc acc sjid workdir; do
    echo "$proc" | tee -a "$accession.log"
    echo "$workdir" | tee -a "$accession.log"
    seff "$sjid" | tee -a "$accession.log"
    cat "$workdir/.command.log" | tee -a "$accession.log"
    echo | tee -a "$accession.log"
  done
echo | tee -a "$accession.log"
done
