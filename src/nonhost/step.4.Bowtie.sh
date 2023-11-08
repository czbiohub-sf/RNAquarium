#!/bin/bash

#SBATCH --job-name=sra_bowtie
#SBATCH --time=14-00:00:00
#SBATCH --array=1-10%10
#SBATCH --partition cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
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

declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

module load anaconda
conda activate "${ENVNAME}"

#remove STAR_out after Bowtie2
remove_STAR_OUT=0

#setting directories
fdir=${WORKING_DIR}/fastq
sdir=${WORKING_DIR}/STAR_out
bdir=${WORKING_DIR}/Bowtie2_out

# declare arrays
readarray -t ACCESSIONS < <(cat "${ACCESSIONS_LIST}") 

bbmap_dir="${TOOLS}/bbmap"
samtools="${TOOLS}/samtools-1.6/samtools"

#setting up index
ZF_idx="${TOOLS}/bowtie2_index/Danio_rerio.GRCz11.dna_sm.primary_assembly"

##########################################################
# check if Bowtie2 result exists, skip current job if so #
##########################################################
# # PE
# if [ -e ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && [ -e ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.gz ] && [ -e ${bdir}/${ACCESSIONS[$idx]}/bowtie2.stats.txt ]
# then 
#     if $(gzip -t ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz) && $(gzip -t ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.gz)
#     then
#         countFileSize=$(gzip -c ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz | wc -c | awk '{print $1}')
#         if [ $countFileSize -gt 100 ]
#         then
#             echo "${ACCESSIONS[$idx]} Bowtie2_out PE result exists and size > 100" >> ${WORKING_DIR}/logs/Bowtie2.skip.log
#             exit 0
#         fi
#     fi
# fi
# # SE
# if [ -e ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && [ -e ${bdir}/${ACCESSIONS[$idx]}/bowtie2.stats.txt ]
# then 
#     if $(gzip -t ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz) 
#     then
#         countFileSize=$(gzip -c ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz | wc -c | awk '{print $1}')
#         if [ $countFileSize -gt 100 ]
#         then
#             echo "${ACCESSIONS[$idx]} Bowtie2_out SE result exists and size > 100" >> ${WORKING_DIR}/logs/Bowtie2.skip.log
#             exit 0
#         fi
#     fi
# fi

###########################
# start executing Bowtie2 #
###########################

if [ -d ${sdir}/PE/${ACCESSIONS[$idx]} ]
then  #paired end
    echo "running bowtie2 with ${ACCESSIONS[$idx]} PE STAR_out" >> ${WORKING_DIR}/logs/Bowtie2.PE.log
    #make bowtie output dir
    rm -rf ${bdir}/${ACCESSIONS[$idx]}
    mkdir ${bdir}/${ACCESSIONS[$idx]}
    #map mate 1 (in unpaired mode)
    if [ -e ${sdir}/PE/${ACCESSIONS[$idx]}/Unmapped.out.mate1.gz ]
    then 
        echo "running bowtie2 with ${ACCESSIONS[$idx]} PE STAR_out Unmapped.out.mate1.gz"
        bowtie2 --quiet --very-sensitive-local \
        --rg-id na --rg LB:na --rg SM:na --rg PL:na --rg PU:na \
        -x ${ZF_idx} -U ${sdir}/PE/${ACCESSIONS[$idx]}/Unmapped.out.mate1.gz | \
        $samtools view -b - >> ${bdir}/${ACCESSIONS[$idx]}/bowtie2.bam
    fi
    #map mate 2 (in unpaired mode)
    if [ -e ${sdir}/PE/${ACCESSIONS[$idx]}/Unmapped.out.mate2.gz ]
    then 
        echo "running bowtie2 with ${ACCESSIONS[$idx]} PE STAR_out Unmapped.out.mate2.gz"
        bowtie2 --quiet --very-sensitive-local \
        --rg-id na --rg LB:na --rg SM:na --rg PL:na --rg PU:na \
        -x ${ZF_idx} -U ${sdir}/PE/${ACCESSIONS[$idx]}/Unmapped.out.mate2.gz | \
        $samtools view -b - >> ${bdir}/${ACCESSIONS[$idx]}/bowtie2.bam
    fi
    #process bam file
    if [ -e ${bdir}/${ACCESSIONS[$idx]}/bowtie2.bam ]
    then #get stats
        echo "processing bowtie2 bam file"
        $samtools view ${bdir}/${ACCESSIONS[$idx]}/bowtie2.bam | cut -f2 | sort | uniq -c > ${bdir}/${ACCESSIONS[$idx]}/bowtie2.stats.txt
        # retrain unmapped reads if bitwise flag contains 4 or 8
        # 4: unmapped 
        # 8: mate unmapped
        $samtools view ${bdir}/${ACCESSIONS[$idx]}/bowtie2.bam | awk '{if( and($2,4)==0 && and($2,8)==0) {print $1}}' > ${bdir}/${ACCESSIONS[$idx]}/bowtie2.mapped.names.txt # print out read name if both read and mate are mapped ( that is SAM bitwise flag does not contain 4 or 8 )
        #filter STAR out using names of mapped reads
        gunzip -c ${sdir}/PE/${ACCESSIONS[$idx]}/Unmapped.out.mate1.gz | ${bbmap_dir}/filterbyname.sh in=stdin names=${bdir}/${ACCESSIONS[$idx]}/bowtie2.mapped.names.txt out=${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT include=f
        gunzip -c ${sdir}/PE/${ACCESSIONS[$idx]}/Unmapped.out.mate2.gz | ${bbmap_dir}/filterbyname.sh in=stdin names=${bdir}/${ACCESSIONS[$idx]}/bowtie2.mapped.names.txt out=${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT include=f
        #compress filtered reads
        gzip ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT
        gzip ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT

        #remove STAR out
        if [ $remove_STAR_OUT -eq 1 ]
        then
            if [ -e ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && $(gzip -t ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz)
            then
                rm ${sdir}/PE/${ACCESSIONS[$idx]}/Unmapped.out.mate1.gz
            fi
            if [ -e ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.gz ] && $(gzip -t ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.gz)
            then
                rm ${sdir}/PE/${ACCESSIONS[$idx]}/Unmapped.out.mate2.gz
            fi
        fi
    #remove bowtie bam
    rm ${bdir}/${ACCESSIONS[$idx]}/bowtie2.bam

    fi
elif [ -d ${sdir}/SE/${ACCESSIONS[$idx]} ]
    then #single end
        echo "running bowtie2 with ${ACCESSIONS[$idx]} SE STAR_out found" >> ${WORKING_DIR}/logs/Bowtie2.SE.log
        #make bowtie output dir
        rm -rf ${bdir}/${ACCESSIONS[$idx]}
        mkdir ${bdir}/${ACCESSIONS[$idx]}

        #map in unpaired mode
        if [ -e ${sdir}/SE/${ACCESSIONS[$idx]}/Unmapped.out.mate1.gz ]
        then 
            echo "running bowtie2 with ${ACCESSIONS[$idx]} SE STAR_out Unmapped.out.mate1.gz"
            bowtie2 --quiet --very-sensitive-local \
            --rg-id na --rg LB:na --rg SM:na --rg PL:na --rg PU:na \
            -x ${ZF_idx} -U ${sdir}/SE/${ACCESSIONS[$idx]}/Unmapped.out.mate1.gz | \
            $samtools view -b - >> ${bdir}/${ACCESSIONS[$idx]}/bowtie2.bam
        fi
        
        #process bam file
        if [ -e ${bdir}/${ACCESSIONS[$idx]}/bowtie2.bam ]
        then #get stats
            echo "processing bowtie2 bam file"
            $samtools view ${bdir}/${ACCESSIONS[$idx]}/bowtie2.bam | cut -f2 | sort | uniq -c > ${bdir}/${ACCESSIONS[$idx]}/bowtie2.stats.txt
            # retain unmapped reads if bitwise flag contains 4 or 8
            # 4: unmapped 
            # 8: mate unmapped
            $samtools view ${bdir}/${ACCESSIONS[$idx]}/bowtie2.bam | awk '{if( and($2,4)==0 && and($2,8)==0) {print $1}}' > ${bdir}/${ACCESSIONS[$idx]}/bowtie2.mapped.names.txt # print out read name if both read and mate are mapped ( that is SAM bitwise flag does not contain 4 or 8 )
            #filter STAR out using names of mapped reads
            gunzip -c ${sdir}/SE/${ACCESSIONS[$idx]}/Unmapped.out.mate1.gz | ${bbmap_dir}/filterbyname.sh in=stdin names=${bdir}/${ACCESSIONS[$idx]}/bowtie2.mapped.names.txt out=${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT include=f
            
            #compress filtered reads
            gzip ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT

            #remove STAR out
            if [ $remove_STAR_OUT -eq 1 ]
            then
                if [ -e ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && $(gzip -t ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz)
                then
                    rm ${sdir}/SE/${ACCESSIONS[$idx]}/Unmapped.out.mate1.gz
                fi
            fi
        fi
        #remove bowtie bam
        rm ${bdir}/${ACCESSIONS[$idx]}/bowtie2.bam

else
    echo "${ACCESSIONS[$idx]} was not found in STAR_out" >> ${WORKING_DIR}/logs/Bowtie2.fail.log
    echo "${ACCESSIONS[$idx]} was not found in STAR_out"
fi