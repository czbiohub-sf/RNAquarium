params.skipHisat = false

/*process parse_fastp {
	cache = false
	input:
	tuple val(idx), val(meta), file("fastp_stats.txt")

	output:
	
	
	script:
	'''
	# FASTP
	fastp_before=\$(sed -n '/Read1 before filtering:/{;n;p;}' fastp_stats.txt | cut -f3 -d' ')
	fastp_after=\$(sed -n '/Read1 after filtering:/{;n;p;}' fastp_stats.txt | cut -f3 -d' ')
	fastp_short=\$(sed -n '/reads failed due to too short:/{;p;}' fastp_stats.txt | cut -f7 -d' ')
	fastp_trimmed=\$(sed -n '/reads with adapter trimmed:/{;p;}' fastp_stats.txt | cut -f5 -d' ')
	printf "%s,%s,%s,%s," \$fastp_before \$fastp_after \$fastp_short \$fastp_trimmed > fastp_stats.csv
	'''
}

process parse_price {
	input:

	output:
	// add env to mates output ?
	tuple val(meta), file("price_stats.csv"), name: csv
	
	script:
	'''
	# PRICE
	# price outputs \b to non-interactive output
	price_total=\$(echo "${price_stats.split("\n")[3].split("/")[1]}" | grep -o "[0-9]\\+")
	price_removed=\$(echo "${price_stats.split("\n")[3].split("/")[2]}" | grep -o "[0-9]\\+")
	price_after=\$(bc <<< "\$price_total - \$price_removed")
	printf "%d,%d," "\$price_total" "\$price_after" > price_stats.csv
	'''
}

process parse_star {

	script:
	'''
	star_before=\$(grep "Number of input reads |" star_stats.txt | grep -o "[0-9]\\+")
	star_avg=\$(grep "Average input read length |" star_stats.txt | grep -o "[0-9]\\+")
	star_unique=\$(grep "Uniquely mapped reads number |" star_stats.txt | grep -o "[0-9]\\+")
	star_multialign=\$(grep "Number of reads mapped to multiple loci |" star_stats.txt | grep -o "[0-9]\\+")
	star_unaligned=\$(bc <<< "\$star_before - \$star_unique - \$star_multialign")
	star_short=\$(grep "Number of reads unmapped: too short |" star_stats.txt | grep -o "[0-9]\\+")
	printf "%s,%s,%s,%s,%s,%s," \$star_before \$star_avg \$star_unique \$star_multialign \$star_unaligned \$star_short
	'''
}*/


process stats_csv {
	cache = false

	input:
	tuple val(idx), val(meta), file("fastp_stats.txt"), val(price_stats), file("hisat2_stats.txt"), file("star_stats.txt"), file("bowtie2_stats.txt"), val(dedup_stats), file(gsnap_stats), val(gsnap_used)

	output:
	stdout

	script:
	def HISAT_UNPAIRED = '''hisat_before=\$(sed -n '/reads; of these:/{;p;}' hisat2_stats.txt | cut -f1 -d' ')
	hisat_unaligned=\$(sed -n '/aligned 0 times/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_aligned_unique=\$(sed -n '/aligned exactly 1 time/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_multialign=\$(sed -n '/aligned >1 times/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_discordant=''
	printf "%s,%s,%s,%s,%s," "\$hisat_before" "\$hisat_unaligned" "\$hisat_aligned_unique" "\$hisat_multialign" "\$hisat_discordant"
	'''
	def HISAT_PAIRED = '''hisat_before=\$(sed -n '/reads; of these:/{;p;}' hisat2_stats.txt | cut -f1 -d' ')
	hisat_unaligned=\$(sed -n '/pairs aligned concordantly 0 times/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_aligned_unique=\$(sed -n '/aligned concordantly exactly 1 time/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_multialign=\$(sed -n '/aligned concordantly >1 times/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	hisat_discordant=\$(sed -n '/aligned discordantly 1 time/{;p;}' hisat2_stats.txt | cut -f5 -d' ')
	printf "%s,%s,%s,%s,%s," "\$hisat_before" "\$hisat_unaligned" "\$hisat_aligned_unique" "\$hisat_multialign" "\$hisat_discordant"
	'''
	def PARSE_HISAT = params.skipHisat ? ",,,,," : (meta.single_end ? HISAT_UNPAIRED : HISAT_PAIRED)
	"""
	#       id single_end     reads       readlen
	printf "id,single_end,starting_reads,median_len,"
	#
	printf "fastp_reads_before,fastp_reads_after,fastp_reads_too_short,fastp_reads_trimmed,"
	#                               %1 - %2
	printf "price_reads_before,price_reads_after,"
	printf "hisat2_reads_before,hisat2_unaligned,hisat2_aligned_unique,hisat2_multialign,hisat2_discordant,"
	printf "star_reads_before,star_avg_len,star_aligned_unique,star_multialign,star_unaligned,star_too_short,"
	printf "bowtie2_reads_before,bowtie2_aligned,bowtie2_multialign,bowtie2_aligned_unique,bowtie2_unaligned,bowtie2_mixed,"
	printf "dedup_reads_before,dedup_reads_after,"
	printf "gsnap_reads_before,gsnap_aligned,gsnap_multialign,gsnap_aligned_unique,gsnap_unaligned,gsnap_mixed,final_reads\n"

	printf "${idx},${meta.single_end},${meta.reads},${meta.readlen},"


	# FASTP
	fastp_before=\$(sed -n '/Read1 before filtering:/{;n;p;}' fastp_stats.txt | cut -f3 -d' ')
	fastp_after=\$(sed -n '/Read1 after filtering:/{;n;p;}' fastp_stats.txt | cut -f3 -d' ')
	fastp_short=\$(sed -n '/reads failed due to too short:/{;p;}' fastp_stats.txt | cut -f7 -d' ')
	fastp_trimmed=\$(sed -n '/reads with adapter trimmed:/{;p;}' fastp_stats.txt | cut -f5 -d' ')
	printf "%s,%s,%s,%s," \$fastp_before \$fastp_after \$fastp_short \$fastp_trimmed

	# PRICE
	# price outputs \b to non-interactive output
	price_total=\$(echo "${price_stats.split("\n")[3].split("/")[1]}" | grep -o "[0-9]\\+")
	price_removed=\$(echo "${price_stats.split("\n")[3].split("/")[2]}" | grep -o "[0-9]\\+")
	price_after=\$(bc <<< "\$price_total - \$price_removed")
	printf "%d,%d," "\$price_total" "\$price_after"

	# HISAT2
	${PARSE_HISAT}

	# STAR
	star_before=\$(grep "Number of input reads |" star_stats.txt | grep -o "[0-9]\\+")
	star_avg=\$(grep "Average input read length |" star_stats.txt | grep -o "[0-9]\\+")
	star_unique=\$(grep "Uniquely mapped reads number |" star_stats.txt | grep -o "[0-9]\\+")
	star_multialign=\$(grep "Number of reads mapped to multiple loci |" star_stats.txt | grep -o "[0-9]\\+")
	star_unaligned=\$(bc <<< "\$star_before - \$star_unique - \$star_multialign")
	star_short=\$(grep "Number of reads unmapped: too short |" star_stats.txt | grep -o "[0-9]\\+")
	printf "%s,%s,%s,%s,%s,%s," \$star_before \$star_avg \$star_unique \$star_multialign \$star_unaligned \$star_short

	# BOWTIE2
	bowtie2_total=0
	bowtie2_multi=0
	bowtie2_aligned=0
	bowtie2_unaligned=0
	bowtie2_mixed=0
	while IFS=' ' read -r count sambits
	do
		# ignore supplementary alignments and second-in-pair
		if [[ \$(( \$sambits & 0x80 || \$sambits & 0x800 )) -ne 0 ]]; then
			continue
		fi
		bowtie2_total=\$(( \$bowtie2_total + \$count ))
		# track secondary mappings but don't add to other counts
		if [[ \$(( "\$sambits" & 0x100 )) -ne 0 ]]; then
			bowtie2_multi=\$(( "\$bowtie2_multi" + \$count ))
			continue
		fi
		if [[ \$(( (\$sambits & 0x5)==0x4 || (\$sambits & 0xD)>0x4 )) -ne 0 ]]; then
			bowtie2_unaligned=\$(( \$bowtie2_unaligned + \$count ))
		fi
		if [[ \$(( ! (\$sambits & 0x4 || \$sambits & 0x8) )) -ne 0 ]]; then
			bowtie2_aligned=\$(( \$bowtie2_aligned + \$count ))
		fi
		if [[ \$(( (\$sambits & 0xD)==0x5 || (\$sambits & 0xD)==0x9  )) -ne 0 ]]; then
			bowtie2_mixed=\$(( \$bowtie2_mixed + \$count ))
		fi
	done < bowtie2_stats.txt
	bowtie2_unique=\$(( \$bowtie2_aligned - \$bowtie2_multi ))
	printf "%s,%s,%s,%s,%s,%s," \$bowtie2_total \$bowtie2_aligned \$bowtie2_multi \$bowtie2_unique \$bowtie2_unaligned \$bowtie2_mixed

	# DEDUP
	dedup_before=\$(echo "$dedup_stats" | grep "total reads:"  | grep -o "[0-9]\\+")
	dedup_after=\$(echo "$dedup_stats" | grep "unique reads:"  | grep -o "[0-9]\\+")
	printf "%s,%s," "\$dedup_before" "\$dedup_after"

	# GSNAP
	if [[ "$gsnap_used" == "no" ]]
	then
	printf "N/A,N/A,N/A,N/A,N/A,N/A,%s\n" "\$dedup_after"
	else
	gsnap_total=0
	gsnap_multi=0
	gsnap_aligned=0
	gsnap_unaligned=0
	gsnap_mixed=0
	while IFS=' ' read -r count sambits
	do
		# ignore supplementary alignments and second-in-pair
		if [[ \$(( \$sambits & 0x80 || \$sambits & 0x800 )) -ne 0 ]]; then
			continue
		fi
		gsnap_total=\$(( \$gsnap_total + \$count ))
		# track secondary mappings but don't add to other counts
		if [[ \$(( "\$sambits" & 0x100 )) -ne 0 ]]; then
			gsnap_multi=\$(( "\$gsnap_multi" + \$count ))
			continue
		fi
		if [[ \$(( (\$sambits & 0x5)==0x4 || (\$sambits & 0xD)>0x4 )) -ne 0 ]]; then
			gsnap_unaligned=\$(( \$gsnap_unaligned + \$count ))
		fi
		if [[ \$(( ! (\$sambits & 0x4 || \$sambits & 0x8) )) -ne 0 ]]; then
			gsnap_aligned=\$(( \$gsnap_aligned + \$count ))
		fi
		if [[ \$(( (\$sambits & 0xD)==0x5 || (\$sambits & 0xD)==0x9  )) -ne 0 ]]; then
			gsnap_mixed=\$(( \$gsnap_mixed + \$count ))
		fi
	done < $gsnap_stats
	gsnap_unique=\$(( \$gsnap_aligned - \$gsnap_multi ))
	printf "%s,%s,%s,%s,%s,%s,%s\n" \$gsnap_total \$gsnap_aligned \$gsnap_multi \$gsnap_unique \$gsnap_unaligned \$gsnap_mixed \$gsnap_unaligned
	fi
	"""
}
