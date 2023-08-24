#!/bin/bash

#SBATCH --job-name=QCfilter
#SBATCH --time=14-00:00:00
#SBATCH --array=1-10%10
#SBATCH --partition cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

#setting directories
WORKING_DIR="${1}"
ACCESSIONS_LIST="${2}"
TOOLS="${3}"
ENVNAME="${4}"
THREAD=4

declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

module load anaconda
conda activate "${ENVNAME}"

# declare arrays
readarray -t ACCESSIONS < <(cat "${ACCESSIONS_LIST}") 

#output directories and paths
FDIR=${WORKING_DIR}/fastq
STDDIR=${WORKING_DIR}/stdout_stderr

fastp_cmd="${TOOLS}/fastp --disable_quality_filtering --disable_length_filtering --compression 6 --thread ${THREAD}"
PriceSeqFilter_cmd="${TOOLS}/PriceSource140408/PriceSeqFilter -a ${THREAD} -rnf 90 -rqf 85 0.98 -log c " # -rqf 85 0.98 was taken from https://github.com/chanzuckerberg/czid-dag/blob/master/idseq_dag/steps/run_priceseq.py
py_file="${tools}/median.py"


############################
# proceed with current job #
############################
echo "QC'ing accession: ${ACCESSIONS[$idx]}" >> ${WORKING_DIR}/logs/QC.process.log
cd $WORKING_DIR

# predict the file names
fq1="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_1.fastq"
fq2="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_2.fastq"
fq="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}.fastq"
trimmed1="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_1.trimmed.fastq"
trimmed2="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_2.trimmed.fastq"
trimmed="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}.trimmed.fastq"
PRICEfiltered1="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_1.PRICEfiltered.fastq"
PRICEfiltered2="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_2.PRICEfiltered.fastq"
PRICEfiltered="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}.PRICEfiltered.fastq"

############################
#check if is scRNAseq data #
############################

if [ -e "${fq1}.gz" ] && [ -e "${fq2}.gz" ]
then
    gunzip ${fq1}.gz ${fq2}.gz
    # use median.py to iterate through fastq file to get median readlength
    length_Read1=$(python ${py_file} ${fq1})
    length_Read2=$(python ${py_file} ${fq1})
    echo "${length_Read1} ${length_Read2}" > "${STDDIR}/${ACCESSIONS[$idx]}/readlength.txt"
    # length_Read1=$(head -n 2 ${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_1.fastq | tail -1 | wc -c)
    # length_Read2=$(head -n 2 ${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_2.fastq | tail -1 | wc -c)
    if [ $length_Read1 -lt 32 ] && [ $length_Read2 -gt 80 ] # read 1 is 26bp cell barcode, read 2 insert size is 90bp
    then
        # echo "${ACCESSIONS[$idx]} scRNAseq ${length_Read1} ${length_Read2}" >> ${working_dir}/logs/SE_PE.log
        echo "${ACCESSIONS[$idx]} scRNAseq" >> "${STDDIR}/${ACCESSIONS[$idx]}/readlength.txt"
        rm ${fq1} #remove read1 (cell barcode)
        mv ${fq2} ${fq}
    else
        # echo "${ACCESSIONS[$idx]} PE ${length_Read1} ${length_Read2}" >> ${working_dir}/logs/SE_PE.log
        echo "${ACCESSIONS[$idx]} PE" >> "${STDDIR}/${ACCESSIONS[$idx]}/readlength.txt"
    fi
else
    if [ -e ${fq}.gz ]
    then
        gunzip ${fq}.gz
        length_Read1=$(head -n 2 ${fq} | tail -1 | wc -c)
        # echo "${ACCESSIONS[$idx]} SE ${length_Read1}" >> ${working_dir}/logs/SE_PE.log
        echo "${ACCESSIONS[$idx]} SE" >> "${STDDIR}/${ACCESSIONS[$idx]}/readlength.txt"
    fi
fi

#########
# fastp #
#########

if [ -e ${fq1} ] && [ -e ${fq2} ]
then
    echo "running fastp in PE mode"
    in1="${fq1}"
    in2="${fq2}"
    out1="${trimmed1}"
    out2="${trimmed2}"
    fastp_cmd+=" --detect_adapter_for_pe --in1 ${in1} --in2 ${in2} --out1 ${out1} --out2 ${out2} --json ${FDIR}/${ACCESSIONS[$idx]}/fastp.json --html ${FDIR}/${ACCESSIONS[$idx]}/fastp.html"
    #run
    $fastp_cmd 1> ${STDDIR}/${ACCESSIONS[$idx]}/fastp.stdout.txt 2> ${STDDIR}/${ACCESSIONS[$idx]}/fastp.stderr.txt

    #remove original fastq dumps
    rm ${fq1}
    rm ${fq2}
fi

if [ -e ${fq} ]
then
    echo "running fastp in SE mode"
    infile="${fq}"
    outfile="${trimmed}"
    fastp_cmd+=" --in1 ${infile} --out1 ${outfile} --json ${FDIR}/${ACCESSIONS[$idx]}/fastp.json --html ${FDIR}/${ACCESSIONS[$idx]}/fastp.html"
    #run
    $fastp_cmd 1> ${STDDIR}/${ACCESSIONS[$idx]}/fastp.stdout.txt 2> ${STDDIR}/${ACCESSIONS[$idx]}/fastp.stderr.txt
    #remove original fastq dumps
    rm ${fq}
fi



##################
# priceseqfilter #
##################
# -a thread, -log c concise output, -fp input, -op output
#-rnf 90:    90 percentage of nucleotides in a read that must be called
#-rqf 85 0.98:       85% of the read length must be 98% accurate (parameterized after CZID, Joe uses 95% of the read length 98% accurate)
if [ -e ${trimmed1} ] && [ -e ${trimmed2} ]
then
    echo "running PriceSeqFilter_cmd"
    PriceSeqFilter_cmd+="-fp ${trimmed1} ${trimmed2} -op ${PRICEfiltered1} ${PRICEfiltered2}"
    #run priceseqfilter
    $PriceSeqFilter_cmd 1> ${STDDIR}/${ACCESSIONS[$idx]}/priceseqfilter.stdout.txt 2> ${STDDIR}/${ACCESSIONS[$idx]}/priceseqfilter.stderr.txt
    #remove trimed fastq
    rm ${trimmed1}
    rm ${trimmed2}
    #rm ${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_1.trimmed.U.fastq
    #rm ${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_2.trimmed.U.fastq
elif [ -e ${trimmed} ]
then
    echo "running PriceSeqFilter_cmd"
    PriceSeqFilter_cmd+="-f ${trimmed} -o ${PRICEfiltered}"
    #run priceseqfilter
    $PriceSeqFilter_cmd 1> ${STDDIR}/${ACCESSIONS[$idx]}/priceseqfilter.stdout.txt 2> ${STDDIR}/${ACCESSIONS[$idx]}/priceseqfilter.stderr.txt
    #remove trimed fastq
    rm ${trimmed}
fi


#compress
if [ -e ${PRICEfiltered1} ]
then
    gzip ${PRICEfiltered1}
fi
if [ -e ${PRICEfiltered2} ]
then
    gzip ${PRICEfiltered2}
fi
if [ -e ${PRICEfiltered} ]
then
    gzip ${PRICEfiltered}
fi
if [ -e ${STDDIR}/${ACCESSIONS[$idx]}/priceseqfilter.stdout.txt ]
then
    gzip ${STDDIR}/${ACCESSIONS[$idx]}/priceseqfilter.stdout.txt
fi
