#!/bin/sh

if [[ $# -lt 1 ]]; then
	echo "usage: collect-reruns.sh input_accessions.txt nonhost_reads_dir .nextflow.log"
	exit 1
fi

# input as either txt or RunInfo csv
if [[ "$1" =~ ".txt" ]]; then
	inlist="$(cat $1)"
else
	# assume RunInfo csv, get first column
	inlist="$(cut -f1 -d',' $1 | tail -n+2)"
fi

# collect input, subtract completed, subtract too short dropouts
outdir=${2:-"nonhost_reads"}
completed="$(ls -1 $outdir)"
nf_log=${3:-".nextflow.log"}
dedup_short=$(grep -o "Process \`dedup (.RR[0-9]\+)\` terminated with an error exit status (1)" $nf_log\
		   | grep -o ".RR[0-9]\+")

# output as txt/csv
if [[ "$1" =~ ".txt" ]]; then
	# plaintext list
	printf '%s\n' $inlist $dedup_short $completed | sort | uniq -u
else
	selected=$(printf '%s\n' $inlist $dedup_short $completed | sort | uniq -u | paste -d"\|")
	head -n1 $1
	grep "$selected" $1
fi

