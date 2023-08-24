#!/bin/bash

#SBATCH --job-name=STAR
#SBATCH --time=14-00:00:00
#SBATCH --array=1-10%10
#SBATCH --partition cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

#setting directories
WORKING_DIR="${1}"
ACCESSIONS_LIST="${2}"
TOOLS="${3}"
ENVNAME="${4}"
REMOVEFQ="${5}"
THREAD=4
count_gene_reads=1 # 1 for run STAR to get read counts for host genes

declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

module load anaconda
conda activate "${ENVNAME}"

# declare arrays
readarray -t ACCESSIONS < <(cat "${ACCESSIONS_LIST}") 

#output directories and paths
FDIR=${WORKING_DIR}/fastq
SDIR=${WORKING_DIR}/STAR_out
STDDIR=${WORKING_DIR}/stdout_stderr
SAMTOOLS_BIN="$TOOLS/samtools-1.16.1/samtools"

echo "running STAR on accession: ${ACCESSIONS[$idx]}" >> ${WORKING_DIR}/logs/STAR.process.log

# predict the file names
PRICEfiltered1gz="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_1.PRICEfiltered.fastq.gz"
PRICEfiltered2gz="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}_2.PRICEfiltered.fastq.gz"
PRICEfilteredgz="${FDIR}/${ACCESSIONS[$idx]}/${ACCESSIONS[$idx]}.PRICEfiltered.fastq.gz"
UMAPPEDMATE1="${SDIR}/PE/${ACCESSIONS[$idx]}/Unmapped.out.mate1"
UMAPPEDMATE2="${SDIR}/PE/${ACCESSIONS[$idx]}/Unmapped.out.mate2"

#generate star command (unmapped)
STAR_bin="${TOOLS}/STAR/STAR-2.7.10a/bin/Linux_x86_64_static/STAR"
indexes_dir="${TOOLS}/STAR/Danio_rerio.GRCz11.108.ERCC" #STAR indexes
indexes_dir2="${TOOLS}/STAR/Danio_rerio.GRCz11.108" #STAR indexes
gtf_noERCC="${TOOLS}/STAR/Danio_rerio.GRCz11.108.gtf"

cmd="${STAR_bin} ";
cmd+="--outFilterMultimapNmax 99999 ";                                       # Joe: capture multimappers
cmd+="--outFilterScoreMinOverLread 0.5 ";                                    # Joe: alignment will be output only if its score is higher than or equal to this value.
cmd+="--outFilterMatchNminOverLread 0.5 ";                                   # Joe: alignment will be output only if the number of matched bases is higher than or equal to this value.
cmd+="--runThreadN ${THREAD} ";                                                      # Duo
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
cmd2+="--runThreadN ${THREAD} ";                                                      # Duo , reduce memory usage, b/c optimizing for gene counts uses more memory
cmd2+="--genomeDir ${indexes_dir2} ";                                         # No ERCC
cmd2+="--readFilesCommand gunzip -c ";                                        # Duo
cmd2+="--limitOutSJcollapsed 200000000 "                                      # Duo: Need to uplift this parameter or else STAR will fail with some runs


#start STAR (paired-end PE), skip if single-ended
if [ -e ${PRICEfiltered1gz} ] && [ -e ${PRICEfiltered2gz} ]
then
    #make STAR output dir
    if [ -e "${SDIR}/PE/${ACCESSIONS[$idx]}" ] 
    then
        rm -rf ${SDIR}/PE/${ACCESSIONS[$idx]}
    fi
    mkdir -p ${SDIR}/PE/${ACCESSIONS[$idx]}

    cmd_PE="${cmd}--outFileNamePrefix ${SDIR}/PE/${ACCESSIONS[$idx]}/ "
    cmd_PE+="--readFilesIn ${PRICEfiltered1gz} ${PRICEfiltered2gz} ";

    #run STAR (unmapped)
    $cmd_PE 1> ${SDIR}/PE/${ACCESSIONS[$idx]}/STAR.stdout.txt 2> ${SDIR}/PE/${ACCESSIONS[$idx]}/STAR.stderr.txt

    #compress 
    if [ -e ${UMAPPEDMATE1} ]
    then
        gzip ${UMAPPEDMATE1}
    fi
    if [ -e ${UMAPPEDMATE2} ]
    then
        gzip ${UMAPPEDMATE2}
    fi
    
    #run STAR (gene read counts)
    if [ ${count_gene_reads} -eq 1 ]
    then 
        #make output dir
        if [ -e "${SDIR}/counts/${ACCESSIONS[$idx]}" ] 
        then
            rm -rf ${SDIR}/counts/${ACCESSIONS[$idx]}
        fi
        mkdir -p ${SDIR}/counts/${ACCESSIONS[$idx]}

        cmd2="${cmd2}--outFileNamePrefix ${SDIR}/counts/${ACCESSIONS[$idx]}/ "
        cmd2+="--readFilesIn ${PRICEfiltered1gz} ${PRICEfiltered2gz} ";
        $cmd2 1> ${SDIR}/counts/${ACCESSIONS[$idx]}/STAR.stdout.txt 2> ${SDIR}/counts/${ACCESSIONS[$idx]}/STAR.stderr.txt
        #Samtools 
        ${SAMTOOLS_BIN} sort -m 45G -n -o ${SDIR}/counts/${ACCESSIONS[$idx]}/Aligned.out.namesorted.bam ${SDIR}/counts/${ACCESSIONS[$idx]}/Aligned.out.bam 2> ${SDIR}/counts/${ACCESSIONS[$idx]}/samtools.stderr.txt
        rm ${SDIR}/counts/${ACCESSIONS[$idx]}/Aligned.out.bam
        #Htseq count
        htseq-count -r name -s no -f bam -m intersection-nonempty ${SDIR}/counts/${ACCESSIONS[$idx]}/Aligned.out.namesorted.bam ${gtf_noERCC} 1> ${SDIR}/counts/${ACCESSIONS[$idx]}/htseq-count.txt 2> ${SDIR}/counts/${ACCESSIONS[$idx]}/htseq-count.stderr.txt
        #remove bam files
        rm ${SDIR}/counts/${ACCESSIONS[$idx]}/Aligned.out.namesorted.bam
    fi

    # remove fastq file
    if [ ${REMOVEFQ} == "yes" ]
    then
        rm -rf ${FDIR}/${ACCESSIONS[$idx]}/*.fastq.gz
        rm -rf ${FDIR}/${ACCESSIONS[$idx]}
    fi
    
    #remove SJ.out.tab
    rm ${SDIR}/PE/${ACCESSIONS[$idx]}/SJ.out.tab
    rm ${SDIR}/PE/${ACCESSIONS[$idx]}/ReadsPerGene.out.tab
    rm ${SDIR}/counts/${ACCESSIONS[$idx]}/SJ.out.tab

##start STAR (single-end and non-paired)
elif [ -e ${PRICEfilteredgz} ] 
then
    #make STAR output dir
    if [ -e "${SDIR}/SE/${ACCESSIONS[$idx]}" ] 
    then
        rm -rf ${SDIR}/SE/${ACCESSIONS[$idx]}
    fi
    mkdir -p ${SDIR}/SE/${ACCESSIONS[$idx]}

    cmd_SE="${cmd}--outFileNamePrefix ${SDIR}/SE/${ACCESSIONS[$idx]}/ "
    cmd_SE+="--readFilesIn ${PRICEfilteredgz} ";

    #run STAR
    $cmd_SE 1> ${SDIR}/SE/${ACCESSIONS[$idx]}/STAR.stdout.txt 2> ${SDIR}/SE/${ACCESSIONS[$idx]}/STAR.stderr.txt

    #compress 
    if [ -e ${SDIR}/SE/${ACCESSIONS[$idx]}/Unmapped.out.mate1 ]
    then
        gzip ${SDIR}/SE/${ACCESSIONS[$idx]}/Unmapped.out.mate1
    fi

    #run STAR (gene read counts)
    if [ ${count_gene_reads} -eq 1 ]
    then 
        #make output dir
        if [ -e "${SDIR}/counts/${ACCESSIONS[$idx]}" ] 
        then
            rm -rf ${SDIR}/counts/${ACCESSIONS[$idx]}
        fi
        mkdir -p ${SDIR}/counts/${ACCESSIONS[$idx]}

        cmd2="${cmd2}--outFileNamePrefix ${SDIR}/counts/${ACCESSIONS[$idx]}/ "
        cmd2+="--readFilesIn ${PRICEfilteredgz} ";
        $cmd2 1> ${SDIR}/counts/${ACCESSIONS[$idx]}/STAR.stdout.txt 2> ${SDIR}/counts/${ACCESSIONS[$idx]}/STAR.stderr.txt
        #Samtools 
        ${SAMTOOLS_BIN} sort -m 45G -n -o ${SDIR}/counts/${ACCESSIONS[$idx]}/Aligned.out.namesorted.bam ${SDIR}/counts/${ACCESSIONS[$idx]}/Aligned.out.bam
        rm ${SDIR}/counts/${ACCESSIONS[$idx]}/Aligned.out.bam
        #Htseq count
        htseq-count -r name -s no -f bam -m intersection-nonempty ${SDIR}/counts/${ACCESSIONS[$idx]}/Aligned.out.namesorted.bam ${gtf_noERCC} > ${SDIR}/counts/${ACCESSIONS[$idx]}/htseq-count.txt 
        #remove bam files
        rm ${SDIR}/counts/${ACCESSIONS[$idx]}/Aligned.out.namesorted.bam
    fi

    # remove fastq file
    if [ ${REMOVEFQ} == "yes" ]
    then
        rm -rf ${FDIR}/${ACCESSIONS[$idx]}/*.fastq.gz
        rm -rf ${FDIR}/${ACCESSIONS[$idx]}
    fi

    #remove SJ.out.tab
    rm ${SDIR}/SE/${ACCESSIONS[$idx]}/SJ.out.tab
    rm ${SDIR}/SE/${ACCESSIONS[$idx]}/ReadsPerGene.out.tab
    rm ${SDIR}/counts/${ACCESSIONS[$idx]}/SJ.out.tab

fi