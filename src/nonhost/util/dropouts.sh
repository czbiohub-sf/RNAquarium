#!/bin/bash
OUTDIR=${1:-"nonhost_reads/"}
NF_LOG=${2:-".nextflow.log"}
download_fail=0
corrupt_file=0
dedup_fail=0
gsnap_fail=$(grep -c "name: gsnap_skip (.*); status: COMPLETED" "$NF_LOG")
other=0
preempt=0
dropout=0
unsuccessful=0
successful=$(ls -1 "$OUTDIR" | wc -l)

mapfile -t fails <<<$(grep -B3 "Error is ignored" ${NF_LOG} | grep "work-dir")

# work directory
for fail in "${fails[@]}"; do
	workdir=$(echo $fail | cut -d"=" -f3) && log=$(cat "$workdir/.command.log")
	if $(echo $fail | grep -q "download"); then
		if [ $VERBOSE ]; then echo "download_fail" $fail &2>/dev/stderr; fi
		((download_fail++))
	elif $(echo $fail | grep -q "dedup"); then
		if [ $VERBOSE ]; then echo "dedup_fail" $fail &2>/dev/stderr; fi
		((dedup_fail++))
	elif $(echo $log | grep -q "PREEMPTED"); then
		if [ $VERBOSE ]; then echo "preempted" $fail &2>/dev/stderr; fi
		((preempt++))
	elif $(echo $fail | grep -q "star_counts"); then
		if [ $VERBOSE ]; then echo "star_counts_fail" $fail &2>/dev/stderr; fi
		continue # STAR_COUNTS, don't double-count
	elif $(echo $log | grep -q ".gz: unexpected end of file|fewer reads in file|more read characters than quality values"); then
		if [ $VERBOSE ]; then echo "corrupte_file" $fail &2>/dev/stderr; fi
		((corrupt_file++))
	else
		if [ $VERBOSE ]; then echo "other" $fail &2>/dev/stderr; fi
		((other++))
	fi
done

unsuccessful=$((download_fail + corrupt_file + preempt + other))
dropout=$((unsuccessful + dedup_fail))
echo "starting_runs,successful,download_fail,corrupt_file,dedup_fail,gsnap_fail,other,preempted,dropout,unsuccessful"
echo "NA,$successful,$download_fail,$corrupt_file,$dedup_fail,$gsnap_fail,$other,$preempt,$dropout,$unsuccessful"
