#!/bin/bash

#SBATCH --job-name=sra_filter_STAR
#SBATCH --time=14-00:00:00
#SBATCH --array=1-54189%100
#SBATCH --partition preempted
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

module load anaconda
conda activate zf_pipeline

#Threads
Trimmomatic_thread=4
PRICEseq_thread=4
STAR_thread=4
fastp_thread=4

#turn on/off gene expr counts
count_gene_reads=1

#use /tmp
use_tmp=0

#setting directories
working_dir="/hpc/scratch/group.theory/jparas/zf_pipeline"


# declare arrays
readarray -t accessions < <(cat ${working_dir}/data/SRA_accession_list.1.27.23.txt) #54189 64G 4cpu
#readarray -t accessions < <(cat /hpc/projects/balla_group/sra_experiments/all_zebrafish_RNAseq/unmapped_dev/accesssions_missing_BT2_res.txt) #
#readarray -t accessions < <(cat /hpc/scratch/group.balla/unmapped_pipeline/STARout.missing.txt)

#output directories
fdir=${working_dir}/fastq
pdir=${working_dir}/prefetched
sdir=${working_dir}/STAR_out

#bins
tools="/hpc/projects/theory_ds/internship/jacob.paras/tools"
sra_bin="${tools}/sratoolkit.3.0.0-ubuntu64/bin"
PRICE_bin="${tools}/PriceSource140408/PriceTI"
STAR_bin="${tools}/STAR/STAR-2.7.10a/bin/Linux_x86_64_static/STAR"
samtools_bin="$tools/samtools-1.16.1/samtools"
Trimmomatic_cmd="${tools}/Trimmomatic-0.38/trimmomatic-0.38.jar"
fastp="${tools}/fastp"
py_file="${tools}/median.py"

#STAR indexes
indexes_dir="${tools}/STAR/Danio_rerio.GRCz11.108.ERCC"
indexes_dir2="${tools}/STAR/Danio_rerio.GRCz11.108"
gtf_noERCC="${tools}/STAR/Danio_rerio.GRCz11.108.gtf"


fastp_cmd="${fastp} --disable_quality_filtering --disable_length_filtering --compression 6 --thread ${fastp_thread}"
PriceSeqFilter_cmd="${tools}/PriceSource140408/PriceSeqFilter -a ${PRICEseq_thread} -rnf 90 -rqf 85 0.98 -log c " # -rqf 85 0.98 was taken from https://github.com/chanzuckerberg/czid-dag/blob/master/idseq_dag/steps/run_priceseq.py

#generate star command (unmapped)
cmd="${STAR_bin} ";
cmd+="--outFilterMultimapNmax 99999 ";                                       # Joe: capture multimappers
cmd+="--outFilterScoreMinOverLread 0.5 ";                                    # Joe: alignment will be output only if its score is higher than or equal to this value.
cmd+="--outFilterMatchNminOverLread 0.5 ";                                   # Joe: alignment will be output only if the number of matched bases is higher than or equal to this value.
cmd+="--runThreadN ${STAR_thread} ";                                                      # Duo
cmd+="--genomeDir ${indexes_dir} ";                                          # Duo
cmd+="--readFilesCommand gunzip -c ";                                        # Duo
cmd+="--outReadsUnmapped Fastx ";                                            # Joe: output of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate file(s).
cmd+="--outFilterMismatchNmax 999 ";                                         # Joe: maximum number of mismatches per pair, large number switches off this filter
cmd+="--outSAMmode None ";                                                   # Joe: disable SAM output
cmd+="--quantMode GeneCounts ";                                              # Joe: output alignments translated into transcript coordinates in the Aligned.toTranscriptome.out.bam file (in addition to alignments in genomic coordinates in Aligned.*.sam/bam files)
cmd+="--clip3pNbases 0 ";                                                    # Joe: This parameter is set to a variable $clip, and I change it to 0
cmd+="--genomeLoad NoSharedMemory ";                                            # Joe:
cmd+="--limitOutSJcollapsed 200000000 "                                      # Duo: Need to uplift this parameter or else STAR will fail with some runs

#generate star command (gene counts)
cmd2="${STAR_bin} ";
#cmd2+="--outFilterType BySJout ";                                       # from QCS  #comment out to reduce memory usage #default: Normal;  BySJout:keep only those reads that contain junctions that passed filtering into SJ.out.tab
cmd2+="--outFilterMultimapNmax 20 ";                                    # from QCS
cmd2+="--alignSJoverhangMin 8 ";                                        # from QCS
cmd2+="--alignSJDBoverhangMin 1 ";                                      # from QCS
cmd2+="--outFilterMismatchNmax 999 ";                                   # from QCS
cmd2+="--outFilterMismatchNoverLmax 0.04 ";                            # from QCS
cmd2+="--alignIntronMin 20 ";                                           # from QCS
cmd2+="--alignIntronMax 1000000 ";                                      # from QCS
cmd2+="--alignMatesGapMax 1000000 ";                                    # from QCS
#cmd2+="--outSAMstrandField intronMotif ";                               # from QCS  #comment out to reduce memory usage, only needed for cufflinks/cuffdiff
cmd2+="--outSAMtype BAM Unsorted ";                                     # from QCS
cmd2+="--outSAMattributes NH HI NM MD ";                                # from QCS
cmd2+="--genomeLoad NoSharedMemory ";                                      # from QCS
cmd2+="--outReadsUnmapped None ";                                       # from QCS
cmd2+="--runThreadN ${STAR_thread} ";                                                      # Duo , reduce memory usage, b/c optimizing for gene counts uses more memory
cmd2+="--genomeDir ${indexes_dir2} ";                                         # No ERCC
cmd2+="--readFilesCommand gunzip -c ";                                        # Duo
cmd2+="--limitOutSJcollapsed 200000000 "                                      # Duo: Need to uplift this parameter or else STAR will fail with some runs


echo "accession: ${accessions[$idx]}"


############################
# proceed with current job #
############################ 
echo "processing..." >> ${working_dir}/logs/STAR.process.log
#sleep 0-120 seconds to space out 100 concurrent jobs from doing prefetch at the same time (prevent timing outs) 
sleep $((idx % 120))
cd $working_dir

#make temp directories if using /tmp
if [ $use_tmp -eq 1 ]
then
    fdir=/tmp/fastq
    pdir=/tmp/prefetched
    sdir=/tmp/STAR_out

    if [ ! -d $fdir ]
    then
        mkdir $fdir
    fi
    
    if [ ! -d $pdir ]
    then
        mkdir $pdir
    fi
    #make STAR-specifc output dir
    if [ ! -d $sdir ]
    then
        mkdir $sdir
        mkdir $sdir/PE
        mkdir $sdir/SE
        mkdir $sdir/counts
        mkdir $sdir/other_stdout_stderr
    fi
fi

# make directories for non STAR stdout and stderr.  STAR out directories for this accession will be purged later in the STAR section
if [ -e "${sdir}/other_stdout_stderr/${accessions[$idx]}" ] 
then
    rm -rf ${sdir}/other_stdout_stderr/${accessions[$idx]}
fi
mkdir ${sdir}/other_stdout_stderr/${accessions[$idx]}

###############
# fasterqdump # 
###############
#create folder for fastq
if [ -e ${fdir}/${accessions[$idx]} ]
then
    rm -rf ${fdir}/${accessions[$idx]}
fi
mkdir ${fdir}/${accessions[$idx]}

#prefetch
${sra_bin}/prefetch --max-size 1t --force ALL --output-directory ${pdir} ${accessions[$idx]} 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/prefetch.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/prefetch.stderr.txt

#fasterq-dump
cd ${pdir}
${sra_bin}/fasterq-dump --split-3 --mem 24G --outdir ${fdir}/${accessions[$idx]} ${accessions[$idx]} 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fqdump.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fqdump.stderr.txt  && {
    echo fasterq-dump: no error
} || {
    echo fasterq-dump encounterd error, reverting back to using fastqdump
    ${sra_bin}/fastq-dump --split-3 --disable-multithreading --outdir ${fdir}/${accessions[$idx]} ${accessions[$idx]} 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fqdump.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/fqdump.stderr.txt
}

#remove prefretched data 
cd ..
rm -rf ${pdir}/${accessions[$idx]}

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

###############
# Trimmomatic #
###############
# SE_adapter=/hpc/projects/balla_group/sra_experiments/tools/Trimmomatic-0.38/adapters/SE.fa
# PE_adapter=/hpc/projects/balla_group/sra_experiments/tools/Trimmomatic-0.38/adapters/PE.fa

# if [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.fastq ] && [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.fastq ]
# then
#     echo "running Trimmomatic PE mode"
#     infiles="${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.fastq ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.fastq"
#     outfiles="${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.trimmed.fastq ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.trimmed.U.fastq ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.trimmed.fastq ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.trimmed.U.fastq"
#     Trimmomatic_cmd+=" PE -threads ${Trimmomatic_thread} -phred33 ${infiles} ${outfiles} ILLUMINACLIP:${PE_adapter}:2:30:10:8:true MINLEN:35"
#     #run 
#     $Trimmomatic_cmd 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/trimmomatic.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/trimmomatic.stderr.txt
#     #remove original fastq dumps
#     rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.fastq
#     rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.fastq
# fi

# if [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}.fastq ]
# then
#     echo "running Trimmomatic SE mode"
#     infiles="${fdir}/${accessions[$idx]}/${accessions[$idx]}.fastq"
#     outfiles="${fdir}/${accessions[$idx]}/${accessions[$idx]}.trimmed.fastq "
#     Trimmomatic_cmd+=" SE -threads ${Trimmomatic_thread} -phred33 ${infiles} ${outfiles} ILLUMINACLIP:${SE_adapter}:2:30:10:8:true MINLEN:35"
#     #run
#     $Trimmomatic_cmd 1> ${sdir}/other_stdout_stderr/${accessions[$idx]}/trimmomatic.stdout.txt 2> ${sdir}/other_stdout_stderr/${accessions[$idx]}/trimmomatic.stderr.txt
#     #remove original fastq dumps
#     rm ${fdir}/${accessions[$idx]}/${accessions[$idx]}.fastq
# fi
# Remove Illumina adapters provided in the fasta file. Initially, look for seed matches
# allowing maximally *2* mismatches. These seeds will be extended and clipped if in the case of paired end
# reads a score of *30* is reached, or in the case of single ended reads a
# score of *10*.
# additional parameters: minAdapterLength = 8, keepBothReads = true; these are set to require pairs to be
#    kept even when an adapter read-through occurs and R2 is a direct reverse complement of R1.


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

#########
# STAR  #
#########

#start STAR (paired-end PE), skip if single-ended
if [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.PRICEfiltered.fastq.gz ] && [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.PRICEfiltered.fastq.gz ]
then
    #make STAR output dir
    if [ -e "${sdir}/PE/${accessions[$idx]}" ] 
    then
        rm -rf ${sdir}/PE/${accessions[$idx]}
    fi
    mkdir ${sdir}/PE/${accessions[$idx]}

    cmd_PE="${cmd}--outFileNamePrefix ${sdir}/PE/${accessions[$idx]}/ "
    cmd_PE+="--readFilesIn ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.PRICEfiltered.fastq.gz ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.PRICEfiltered.fastq.gz ";

    #run STAR (unmapped)
    $cmd_PE 1> ${sdir}/PE/${accessions[$idx]}/STAR.stdout.txt 2> ${sdir}/PE/${accessions[$idx]}/STAR.stderr.txt

    #compress 
    if [ -e ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate1 ]
    then
        gzip ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate1
    fi
    if [ -e ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate2 ]
    then
        gzip ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate2
    fi
    
    #run STAR (gene read counts)
    if [ ${count_gene_reads} -eq 1 ]
    then 
        #make output dir
        if [ -e "${sdir}/counts/${accessions[$idx]}" ] 
        then
            rm -rf ${sdir}/counts/${accessions[$idx]}
        fi
        mkdir ${sdir}/counts/${accessions[$idx]}

        cmd2="${cmd2}--outFileNamePrefix ${sdir}/counts/${accessions[$idx]}/ "
        cmd2+="--readFilesIn ${fdir}/${accessions[$idx]}/${accessions[$idx]}_1.PRICEfiltered.fastq.gz ${fdir}/${accessions[$idx]}/${accessions[$idx]}_2.PRICEfiltered.fastq.gz ";
        $cmd2 1> ${sdir}/counts/${accessions[$idx]}/STAR.stdout.txt 2> ${sdir}/counts/${accessions[$idx]}/STAR.stderr.txt
        #Samtools 
        ${samtools_bin} sort -m 45G -n -o ${sdir}/counts/${accessions[$idx]}/Aligned.out.namesorted.bam ${sdir}/counts/${accessions[$idx]}/Aligned.out.bam 2> ${sdir}/counts/${accessions[$idx]}/samtools.stderr.txt
        rm ${sdir}/counts/${accessions[$idx]}/Aligned.out.bam
        #Htseq count
        htseq-count -r name -s no -f bam -m intersection-nonempty ${sdir}/counts/${accessions[$idx]}/Aligned.out.namesorted.bam ${gtf_noERCC} 1> ${sdir}/counts/${accessions[$idx]}/htseq-count.txt 2> ${sdir}/counts/${accessions[$idx]}/htseq-count.stderr.txt
        #remove bam files
        rm ${sdir}/counts/${accessions[$idx]}/Aligned.out.namesorted.bam
    fi

    # remove fastq file
    rm -rf ${fdir}/${accessions[$idx]}/*.fastq.gz
    rm -rf ${fdir}/${accessions[$idx]}
    
    #remove SJ.out.tab
    rm ${sdir}/PE/${accessions[$idx]}/SJ.out.tab
    rm ${sdir}/PE/${accessions[$idx]}/ReadsPerGene.out.tab
    rm ${sdir}/counts/${accessions[$idx]}/SJ.out.tab

##start STAR (single-end and non-paired)
elif [ -e ${fdir}/${accessions[$idx]}/${accessions[$idx]}.PRICEfiltered.fastq.gz ] 
then
    #make STAR output dir
    if [ -e "${sdir}/SE/${accessions[$idx]}" ] 
    then
        rm -rf ${sdir}/SE/${accessions[$idx]}
    fi
    mkdir ${sdir}/SE/${accessions[$idx]}

    cmd_SE="${cmd}--outFileNamePrefix ${sdir}/SE/${accessions[$idx]}/ "
    cmd_SE+="--readFilesIn ${fdir}/${accessions[$idx]}/${accessions[$idx]}.PRICEfiltered.fastq.gz ";

    #run STAR
    $cmd_SE 1> ${sdir}/SE/${accessions[$idx]}/STAR.stdout.txt 2> ${sdir}/SE/${accessions[$idx]}/STAR.stderr.txt

    #compress 
    if [ -e ${sdir}/SE/${accessions[$idx]}/Unmapped.out.mate1 ]
    then
        gzip ${sdir}/SE/${accessions[$idx]}/Unmapped.out.mate1
    fi

    #run STAR (gene read counts)
    if [ ${count_gene_reads} -eq 1 ]
    then 
        #make output dir
        if [ -e "${sdir}/counts/${accessions[$idx]}" ] 
        then
            rm -rf ${sdir}/counts/${accessions[$idx]}
        fi
        mkdir ${sdir}/counts/${accessions[$idx]}

        cmd2="${cmd2}--outFileNamePrefix ${sdir}/counts/${accessions[$idx]}/ "
        cmd2+="--readFilesIn ${fdir}/${accessions[$idx]}/${accessions[$idx]}.PRICEfiltered.fastq.gz ";
        $cmd2 1> ${sdir}/counts/${accessions[$idx]}/STAR.stdout.txt 2> ${sdir}/counts/${accessions[$idx]}/STAR.stderr.txt
        #Samtools 
        ${samtools_bin} sort -m 45G -n -o ${sdir}/counts/${accessions[$idx]}/Aligned.out.namesorted.bam ${sdir}/counts/${accessions[$idx]}/Aligned.out.bam
        rm ${sdir}/counts/${accessions[$idx]}/Aligned.out.bam
        #Htseq count
        htseq-count -r name -s no -f bam -m intersection-nonempty ${sdir}/counts/${accessions[$idx]}/Aligned.out.namesorted.bam ${gtf_noERCC} > ${sdir}/counts/${accessions[$idx]}/htseq-count.txt 
        #remove bam files
        rm ${sdir}/counts/${accessions[$idx]}/Aligned.out.namesorted.bam
    fi

    # remove fastq file
    rm -rf ${fdir}/${accessions[$idx]}/*.fastq.gz
    rm -rf ${fdir}/${accessions[$idx]}

    #remove SJ.out.tab
    rm ${sdir}/SE/${accessions[$idx]}/SJ.out.tab
    rm ${sdir}/SE/${accessions[$idx]}/ReadsPerGene.out.tab
    rm ${sdir}/counts/${accessions[$idx]}/SJ.out.tab

fi

# if using /tmp, copy result files out of /tmp
if [ ${use_tmp} -eq 1 ]
then
    sdir2=${working_dir}/STAR_out

    if [ -e "${sdir}/PE/${accessions[$idx]}" ]
    then
        rm -rf ${sdir2}/PE/${accessions[$idx]}
        mv ${sdir}/PE/${accessions[$idx]} ${sdir2}/PE/${accessions[$idx]}
    fi

    if [ -e "${sdir}/SE/${accessions[$idx]}" ]
    then
        rm -rf ${sdir2}/SE/${accessions[$idx]}
        mv ${sdir}/SE/${accessions[$idx]} ${sdir2}/SE/${accessions[$idx]}
    fi

    if [ -e "${sdir}/counts/${accessions[$idx]}" ]
    then
        rm -rf ${sdir2}/counts/${accessions[$idx]}
        mv ${sdir}/counts/${accessions[$idx]} ${sdir2}/counts/${accessions[$idx]}
    fi

    if [ -e "${sdir}/other_stdout_stderr/${accessions[$idx]}" ]
    then
        rm -rf ${sdir2}/other_stdout_stderr/${accessions[$idx]}
        mv ${sdir}/other_stdout_stderr/${accessions[$idx]} ${sdir2}/other_stdout_stderr/${accessions[$idx]}
    fi

fi
