#/usr/bin/env bash

set +m
shopt -s lastpipe

printf "acc,seq_type,mate2_median\n"
while IFS= read -r dir;
do
	printf "%s,%s,%s\n" $(cat "$dir/.command.out") $(grep median "$dir/.command.env" | cut -d'=' -f2,2)
done
