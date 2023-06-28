############################
#check if is scRNAseq data #
############################
if [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.fastq ] && [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.fastq ]
then
    # use median.py to iterate through fastq file to get median readlength
    length_Read1=$(python ${py_file} ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.fastq)
    length_Read2=$(python ${py_file} ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.fastq)
    echo "${length_Read1} ${length_Read2}" > "${sdir}/other_stdout_stderr/${accessions[$idx]}/readlength.txt"
    # length_Read1=$(head -n 2 ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.fastq | tail -1 | wc -c)
    # length_Read2=$(head -n 2 ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.fastq | tail -1 | wc -c)
    if [ $length_Read1 -lt 32 ] && [ $length_Read2 -gt 80 ] # read 1 is 26bp cell barcode, read 2 insert size is 90bp
    then
        # echo "${accessions[$idx]} scRNAseq ${length_Read1} ${length_Read2}" >> ${working_dir}/logs/SE_PE.log
        echo "${accessions[$idx]} scRNAseq" >> "${sdir}/other_stdout_stderr/${accessions[$idx]}/readlength.txt"
        rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.fastq #remove read1 (cell barcode)
        mv ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.fastq ${fdir}/${accessions[$idx]}/${accessions[$idx]}.fastq
    else
        # echo "${accessions[$idx]} PE ${length_Read1} ${length_Read2}" >> ${working_dir}/logs/SE_PE.log
        echo "${accessions[$idx]} PE" >> "${sdir}/other_stdout_stderr/${accessions[$idx]}/readlength.txt"
    fi
else
    if [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}.fastq ]
    then
        length_Read1=$(head -n 2 ${fdir}/${accessions[$idx]}/${accessions[$idx]}.fastq | tail -1 | wc -c)
        # echo "${accessions[$idx]} SE ${length_Read1}" >> ${working_dir}/logs/SE_PE.log
        echo "${accessions[$idx]} SE" >> "${sdir}/other_stdout_stderr/${accessions[$idx]}/readlength.txt"
    fi
fi

#########
# fastp #
#########

if [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.fastq ] && [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.fastq ]
then
    echo "running fastp in PE mode"
    in1="${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.fastq"
    in2="${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.fastq"
    out1="${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.trimmed.fastq"
    out2="${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.trimmed.fastq"
    fastp_cmd+=" --detect_adapter_for_pe --in1 ${in1} --in2 ${in2} --out1 ${out1} --out2 ${out2} --json ${fdir}/${accessions[$idx]}/fastp.json --html ${fdir}/${accessions[$idx]}/fastp.html"
    #run
    $fastp_cmd 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fastp.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fastp.stderr.txt

    #remove original fastq dumps
    rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.fastq
    rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.fastq
fi

if [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}.fastq ]
then
    echo "running fastp in SE mode"
    infile="${fdir}/${accessions[$idx]}/${accessions[$idx]}.fastq"
    outfile="${fdir}/${accessions[$idx]}/${accessions[$idx]}.trimmed.fastq"
    fastp_cmd+=" --in1 ${infile} --out1 ${outfile} --json ${fdir}/${accessions[$idx]}/fastp.json --html ${fdir}/${accessions[$idx]}/fastp.html"
    #run
    $fastp_cmd 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fastp.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fastp.stderr.txt
    #remove original fastq dumps
    rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}.fastq
fi



##################
# priceseqfilter #
##################
# -a thread, -log c concise output, -fp input, -op output
#-rnf 90:    90 percentage of nucleotides in a read that must be called
#-rqf 85 0.98:       85% of the read length must be 98% accurate (parameterized after CZID, Joe uses 95% of the read length 98% accurate)
if [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.trimmed.fastq ] && [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.trimmed.fastq ]
then
    echo "running PriceSeqFilter_cmd"
    PriceSeqFilter_cmd+="-fp ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.trimmed.fastq ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.trimmed.fastq -op ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.PRICEfiltered.fastq ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.PRICEfiltered.fastq"
    #run priceseqfilter
    $PriceSeqFilter_cmd 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/priceseqfilter.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/priceseqfilter.stderr.txt
    #remove trimed fastq
    rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.trimmed.fastq
    rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.trimmed.fastq
    #rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.trimmed.U.fastq
    #rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.trimmed.U.fastq
elif [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}.trimmed.fastq ]
then
    echo "running PriceSeqFilter_cmd"
    PriceSeqFilter_cmd+="-f ${fdir}/${accessions[$idx]}/${accessions[$idx]}.trimmed.fastq -o ${fdir}/${accessions[$idx]}/${accessions[$idx]}.PRICEfiltered.fastq"
    #run priceseqfilter
    $PriceSeqFilter_cmd 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/priceseqfilter.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/priceseqfilter.stderr.txt
    #remove trimed fastq
    rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}.trimmed.fastq
fi


#compress
if [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.PRICEfiltered.fastq ]
then
    gzip ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.PRICEfiltered.fastq
fi
if [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.PRICEfiltered.fastq ]
then
    gzip ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.PRICEfiltered.fastq
fi
if [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}.PRICEfiltered.fastq ]
then
    gzip ${fdir}/${accessions[$idx]}/${accessions[$idx]}.PRICEfiltered.fastq
fi
if [ -e ${sdir}/other_stdout_stderr/${accessions[$idx]}/priceseqfilter.stdout.txt ]
then
    gzip ${sdir}/other_stdout_stderr/${accessions[$idx]}/priceseqfilter.stdout.txt
fi
