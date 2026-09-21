// ===========================================================================
// WALKTHROUGH-LOCAL PATCH of src/nonhost/modules/local/stats.nf
// Makes the stats row robust to empty/missing per-step values on tiny test
// inputs, so it always has its full 41 fields and the end-of-run summary
// reduce in main.nf never sees a null (which otherwise NPEs). Two failure
// modes are covered: (1) $(( )) / bc subtractions on EMPTY operands (e.g.
// the kallisto kb_unaligned line -> 'operand expected'); (2) UNQUOTED empty
// printf args being dropped by word-splitting -> short row. All arithmetic is
// defaulted to 0 and all printf args are quoted+defaulted; parse_sam also
// skips a missing/empty samfile. NO counts/math change for present, non-empty
// values. Test-drive convenience only; NOT part of the upstream PR.
// (Note: the kallisto branch guard `[[ "input.2" == "na" ]]` compares a
//  FILENAME string, not file content -- a likely upstream bug worth a look.)
// ===========================================================================
params.skipHisat = false

process stats_csv {
	cache = false

	input:
	tuple val(idx), val(meta), file("fastp_stats.txt"), file(kb_stats), file("hisat2_stats.txt"), file("star_stats.txt"), file("bowtie2_stats.txt"), val(dedup_stats), file(gsnap_stats), val(gsnap_used)

	output:
	stdout

	script:
	"""
parse_sam() {
	samfile="\$1"
	final_count=\${2:-}
	total=0
	multi=0
	aligned=0
	unaligned=0
	mixed=0
	unique=0
	if [[ -s "\$samfile" ]]; then
	while IFS=' ' read -r count sambits
	do
	    # ignore supplementary alignments and second-in-pair
	    #  SECONDINPAIR || QCFAIL || SUPPLEMENTARY 
		if [[ \$(( \$sambits & 0x80 || \$sambits & 0x200 || \$sambits & 0x800 )) -ne 0 ]]; then
			continue
		fi
		# track secondary mappings but don't add to other counts
        # SECONDARY
	    if [[ \$(( "\$sambits" & 0x100 )) -ne 0 ]]; then
			multi=\$(( "\$multi" + \$count ))
			continue
	    fi
		total=\$(( \$total + \$count ))
        # [SE] PAIR|UNMAP == UNMAP        || [PE] PAIR|UNMAP|MUNMAP > UNMAP [PAIR+UNMAP,PAIR+MUNMAP,PAIR+UNMAP+MUNMAP]
		if [[ \$(( (\$sambits & 0x5)==0x4 || (\$sambits & 0xD)>0x4 )) -ne 0 ]]; then
			unaligned=\$(( \$unaligned + \$count ))
	    fi
	    # [CONCORDANT] ! (UNMAP || MUNMAP)
    	if [[ \$(( ! ( (\$sambits & 0x4) || (\$sambits & 0x8) ) )) -ne 0 ]]; then
			aligned=\$(( \$aligned + \$count ))
	    fi
	    # PAIR|UNMAP|MUNMAP == PAIR|UNMAP || PAIR|UNMAP|MUNMAP == PAIR|MUNMAP
		if [[ \$(( (\$sambits & 0xD)==0x5 || (\$sambits & 0xD)==0x9  )) -ne 0 ]]; then
			mixed=\$(( \$mixed + \$count ))
		fi
	done < "\$samfile"
	fi
	unique=\$(( \${aligned:-0} - \${multi:-0} ))
	printf "%s,%s,%s,%s,%s,%s," "\${total:-0}" "\${aligned:-0}" "\${multi:-0}" "\${unique:-0}" "\${unaligned:-0}" "\${mixed:-0}"
	if [[ -n "\$final_count" ]]; then
		printf "%s\\n" "\${unaligned:-0}"
	fi
}
	
	#       id single_end     reads       readlen
	printf "id,single_end,starting_reads,r1_median_len,r2_median_len,"
	#
	printf "fastp_reads_before,fastp_reads_after,fastp_reads_too_short,fastp_reads_trimmed,"
	#                                                                           %1 - %2
	printf "kallisto_reads_before,kallisto_aligned,kallisto_aligned_unique,kallisto_unaligned,kallisto_targets,"
	#                               %1 - %2
	printf "hisat2_reads_before,hisat2_aligned,hisat2_multialign,hisat2_aligned_unique,hisat2_unaligned,hisat2_mixed,"
	printf "star_reads_before,star_avg_len,star_aligned_unique,star_multialign,star_unaligned,star_too_short,"
	printf "bowtie2_reads_before,bowtie2_aligned,bowtie2_multialign,bowtie2_aligned_unique,bowtie2_unaligned,bowtie2_mixed,"
	printf "dedup_reads_before,dedup_reads_after,"
	printf "gsnap_reads_before,gsnap_aligned,gsnap_multialign,gsnap_aligned_unique,gsnap_unaligned,gsnap_mixed,final_reads\n"

	printf "${idx},${meta.single_end},${meta.reads},${meta.readlen},${meta.readlen_2},"


	# FASTP
	fastp_before=\$(sed -n '/Read1 before filtering:/{;n;p;}' fastp_stats.txt | cut -f3 -d' ')
	fastp_after=\$(sed -n '/Read1 after filtering:/{;n;p;}' fastp_stats.txt | cut -f3 -d' ')
	fastp_short=\$(sed -n '/reads failed due to too short:/{;p;}' fastp_stats.txt | cut -f7 -d' ')
	fastp_trimmed=\$(sed -n '/reads with adapter trimmed:/{;p;}' fastp_stats.txt | cut -f5 -d' ')
	printf "%s,%s,%s,%s," "\${fastp_before:-0}" "\${fastp_after:-0}" "\${fastp_short:-0}" "\${fastp_trimmed:-0}"

	# KALLISTO
	if [[ "$kb_stats" == "na" ]]
	then
		printf "N/A,N/A,N/A,N/A,N/A,"
	else
		kb_targets=\$(sed -nE 's/\t+"n_targets": +([0-9]+),/\\1/p' "$kb_stats")
		kb_before=\$(sed -nE 's/\t+"n_processed": +([0-9]+),/\\1/p' "$kb_stats")
		kb_aligned=\$(sed -nE 's/\t+"n_pseudoaligned": +([0-9]+),/\\1/p' "$kb_stats")
		kb_unique=\$(sed -nE 's/\t+"n_unique": +([0-9]+),/\\1/p' "$kb_stats")
		kb_unaligned=\$(( \${kb_before:-0} - \${kb_aligned:-0} ))
		printf "%s,%s,%s,%s,%s," "\${kb_before:-0}" "\${kb_aligned:-0}" "\${kb_unique:-0}" "\${kb_unaligned:-0}" "\${kb_targets:-0}"
	fi

	# HISAT2
	parse_sam "hisat2_stats.txt"

	# STAR
	# Log.final.out gives a unique field and is used elsewhere, so not sam_stats
	star_before=\$(grep "Number of input reads |" star_stats.txt | grep -o "[0-9]\\+")
	star_avg=\$(grep "Average input read length |" star_stats.txt | grep -o "[0-9]\\+")
	star_unique=\$(grep "Uniquely mapped reads number |" star_stats.txt | grep -o "[0-9]\\+")
	star_multialign=\$(grep "Number of reads mapped to multiple loci |" star_stats.txt | grep -o "[0-9]\\+")
	star_unaligned=\$(bc <<< "\${star_before:-0} - \${star_unique:-0} - \${star_multialign:-0}")
	star_short=\$(grep "Number of reads unmapped: too short |" star_stats.txt | grep -o "[0-9]\\+")
	printf "%s,%s,%s,%s,%s,%s," "\${star_before:-0}" "\${star_avg:-0}" "\${star_unique:-0}" "\${star_multialign:-0}" "\${star_unaligned:-0}" "\${star_short:-0}"
	
	# BOWTIE2
	parse_sam "bowtie2_stats.txt"
	
	# DEDUP
	dedup_before=\$(echo "$dedup_stats" | grep "total reads:"  | grep -o "[0-9]\\+")
	dedup_after=\$(echo "$dedup_stats" | grep "unique reads:"  | grep -o "[0-9]\\+")
	printf "%s,%s," "\${dedup_before:-0}" "\${dedup_after:-0}"

	# GSNAP
	if [[ "$gsnap_used" == "no" ]]
	then
		printf "N/A,N/A,N/A,N/A,N/A,N/A,%s\\n" "\$dedup_after"
	else
		parse_sam "$gsnap_stats" 1
	fi
	"""
}
