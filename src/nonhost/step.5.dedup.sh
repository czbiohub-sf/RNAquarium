#!/bin/bash

#SBATCH --job-name=sra_dedup
#SBATCH --time=14-00:00:00
#SBATCH --array=1-10%10
#SBATCH --nodes=1
#SBATCH --partition cpu
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
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

#parameters
perc_len=100 #percent length used for dedup
rm_BT2_res=0 #whether to remove bowtie2 results, 0=don't remove, 1=remove
minimum_read_length=20 #discard *runs* with read length shorter than this threshold

#setting directories
bdir=${WORKING_DIR}/Bowtie2_out
ddir=${WORKING_DIR}/Dedup_out
sdir=${WORKING_DIR}/STAR_out

# declare arrays
readarray -t ACCESSIONS < <(cat "${ACCESSIONS_LIST}") 

czid_dedup_bin_path="${TOOLS}/czid-dedup-Linux"

#main command

#check if result already exists for the current accession, if so, check gzip file integrity. If all pass, exit the script
if [ -e ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz ] && [ -e ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq.gz ] #check if dedup result files exist, PE reads
then
    if $(gzip -t ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz) && $(gzip -t ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq.gz) #check if dedup result files pass gzip integrity test, PE reads
        then
            countFileSize1=$(gzip -c ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz | wc -c | awk '{print $1}')
            countFileSize2=$(gzip -c ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq.gz | wc -c | awk '{print $1}')
            if  [ $countFileSize1 -gt 100 ] && [ $countFileSize2 -gt 100 ]
            then
                echo "Skipping ${ddir}/${ACCESSIONS[$idx]} PE reads, dedup.fastq.gz found and size > 100" >> ${WORKING_DIR}/logs/dedup.skip.log
                exit 0
            fi
    fi
else #check SE
    if [ -e ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz ] #check if dedup result files exist, SE
    then
        if $(gzip -t ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz) #check if dedup result files pass gzip integrity test, SE
            then
                countFileSize1=$(gzip -c ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz | wc -c | awk '{print $1}')
                if  [ $countFileSize1 -gt 100 ]
                then
                    echo "skipping ${ddir}/${ACCESSIONS[$idx]} SE reads, dedup.fastq.gz found and size > 100" >> ${WORKING_DIR}/logs/dedup.skip.log
                    exit 0
                fi
        fi
    fi
fi


#########################################
#Dedup results not found or is corrupted#
#Start fresh dedup for current accession#
#########################################
# PE reads
if [ -e ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && [ -e ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.gz ] && [ $(gzip -l ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz | awk 'NR==2 {print $2}') -ne 0 ] && [ $(gzip -l ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.gz | awk 'NR==2 {print $2}') -ne 0 ]  #gzip file exists and is not empty
then 
    if $(gzip -t ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz) && $(gzip -t ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.gz)
    then
        echo "running dedup using ${bdir}/${ACCESSIONS[$idx]} PE reads" >> ${WORKING_DIR}/logs/dedup.PE.log
        if [ ! -d "${ddir}/${ACCESSIONS[$idx]}" ] 
        then
            mkdir ${ddir}/${ACCESSIONS[$idx]}
        fi
        gunzip -c ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz > ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq
        gunzip -c ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.gz > ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.fastq

        #determine read length
        read_length=0
        if [ -e ${sdir}/PE/${ACCESSIONS[$idx]}/Log.final.out ]
        then
            read_length_from_STAR_report=$(cat ${sdir}/PE/${ACCESSIONS[$idx]}/Log.final.out | sed '7q;d'| cut -f2)
            read_length=$read_length_from_STAR_report
        else
            echo "STAR final log not found: ${sdir}/PE/${ACCESSIONS[$idx]}/Log.final.out, determining read length with ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq" >> ${WORKING_DIR}/logs/dedup.warning.log
            read_length_from_fastq=$(sed '2q;d' ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq | wc -m)
            read_length=$read_length_from_fastq
        fi
        re='^[0-9]+$'
        read_length=$(echo $read_length | sed 's/[^0-9]*//g')
        if [ $read_length -eq 0 ] || ! [[ $read_length =~ $re ]] || (( $read_length < $minimum_read_length ))
        then
            echo "[Fatal error][${ACCESSIONS[$idx]}] Cannot determine read length or read too short, got read length: $read_length, " >> ${WORKING_DIR}/logs/dedup.error.log
            exit 1
        fi
        length_useFor_dedup=$(bc <<< "scale=3; $(bc <<< "scale=3; $read_length / 100") * $perc_len") # calculate percent read length
        length_useFor_dedup=${length_useFor_dedup%.*}
        echo "${ACCESSIONS[$idx]} read length ${read_length}, use ${length_useFor_dedup} for dedup">> ${WORKING_DIR}/logs/dedup.lengths.log

        #run dedup with -l percent read length
        ${czid_dedup_bin_path} \
        -i ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq \
        -i ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.fastq \
        -o ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq \
        -o ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq \
        -c ${ddir}/${ACCESSIONS[$idx]}/cluster.csv \
        -l ${length_useFor_dedup}
        rm ${ddir}/${ACCESSIONS[$idx]}/cluster.csv #  writing to cluster.csv then delete it is to preventing multiple jobs trying to write to the same file and causing unforseen problems
        rm ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq #remove unzipped BT2 results
        rm ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.fastq #remove unzipped BT2 results
        gzip ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq
        gzip ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq
        if [ $rm_BT2_res -eq 1 ] #remove bowtie2 results
        then    
            rm ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz
            rm ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate2.filteredbyBT.gz
        fi
    else
        echo "[Fatal error][${ACCESSIONS[$idx]}] At least one of the ${bdir}/${ACCESSIONS[$idx]}/*.filteredbyBT.gz files failed the integrity test" >> ${WORKING_DIR}/logs/dedup.error.log
    fi
else #SE reads
    if [ -e ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && [ $(gzip -l ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz | awk 'NR==2 {print $2}') -ne 0 ] #gzip file exists and is not empty
    then
        if $(gzip -t ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz)
        then
            echo "running dedup using ${bdir}/${ACCESSIONS[$idx]} SE reads" >> ${WORKING_DIR}/logs/dedup.SE.log
            if [ ! -d "${ddir}/${ACCESSIONS[$idx]}" ] 
            then
                mkdir ${ddir}/${ACCESSIONS[$idx]}
            fi
            gunzip -c ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz > ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq

            #determine read length
            read_length=0
            if [ -e ${sdir}/SE/${ACCESSIONS[$idx]}/Log.final.out ]
            then
                read_length_from_STAR_report=$(cat ${sdir}/SE/${ACCESSIONS[$idx]}/Log.final.out | sed '7q;d'| cut -f2)
                read_length=$read_length_from_STAR_report
            else
                echo "STAR final log not found: ${sdir}/SE/${ACCESSIONS[$idx]}/Log.final.out, determining read length with ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq" >> ${WORKING_DIR}/logs/dedup.warning.log
                read_length_from_fastq=$(sed '2q;d' ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq | wc -m)
                read_length=$read_length_from_fastq
            fi
            re='^[0-9]+$'
            read_length=$(echo $read_length | sed 's/[^0-9]*//g')
            if [ $read_length -eq 0 ] || ! [[ $read_length =~ $re ]] || (( $read_length < $minimum_read_length ))
            then
                echo "[Fatal error][${ACCESSIONS[$idx]}] Cannot determine read length or read too short, got read length: $read_length, " >> ${WORKING_DIR}/logs/dedup.error.log
                exit 1
            fi
            length_useFor_dedup=$(bc <<< "scale=3; $(bc <<< "scale=3; $read_length / 100") * $perc_len") # calculate percent read length
            length_useFor_dedup=${length_useFor_dedup%.*}
            echo "${ACCESSIONS[$idx]} read length ${read_length}, use ${length_useFor_dedup} for dedup">> ${WORKING_DIR}/logs/dedup.lengths.log

            #run dedup with -l percent read length
            ${czid_dedup_bin_path} \
            -i ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq \
            -o ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq \
            -c ${ddir}/${ACCESSIONS[$idx]}/cluster.csv \
            -l ${length_useFor_dedup}
            rm ${ddir}/${ACCESSIONS[$idx]}/cluster.csv # writing to cluster.csv then delete it is to preventing multiple jobs trying to write to the same file and causing unforseen problems
            rm ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq #remove unzipped BT2 results
            gzip ${ddir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq
            if [ $rm_BT2_res -eq 1 ] #remove bowtie2 results
            then    rm ${bdir}/${ACCESSIONS[$idx]}/Unmapped.out.mate1.filteredbyBT.gz
            fi
        else
            echo "[Fatal error][${ACCESSIONS[$idx]}] At least one of the ${bdir}/${ACCESSIONS[$idx]}/*.filteredbyBT.gz files failed the integrity test" >> ${WORKING_DIR}/logs/dedup.error.log
        fi
    else 
        echo "[Fatal error][${ACCESSIONS[$idx]}] At least one of the ${bdir}/${ACCESSIONS[$idx]}/*.filteredbyBT.gz files were not found or are empty files" >> ${WORKING_DIR}/logs/dedup.error.log
    fi
fi
###########
#end dedup#
###########




