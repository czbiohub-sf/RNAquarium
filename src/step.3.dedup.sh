#!/bin/bash

#SBATCH --job-name=sra_dedup
#SBATCH --time=14-00:00:00
#SBATCH --array=1-54189%500
#SBATCH --nodes=1
#SBATCH --partition preempted
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH -e slurm.out.dedup/slurm-%A_%a.err
#SBATCH -o slurm.out.dedup/slurm-%A_%a.out

declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

module load anaconda
conda activate zf_pipeline

#parameters
perc_len=100 #percent length used for dedup
rm_BT2_res=0 #whether to remove bowtie2 results, 0=don't remove, 1=remove
minimum_read_length=20 #discard *runs* with read length shorter than this threshold

#setting directories
working_dir="/hpc/scratch/group.theory/jparas/zf_pipeline"
bdir=${working_dir}/Bowtie2_out
ddir=${working_dir}/Dedup_out
sdir=${working_dir}/STAR_out

# declare arrays
readarray -t accessions < <(cat "${working_dir}/data/SRA_accession_list.1.27.23.txt")

tools="/hpc/projects/theory/sharing/internship/jacob.paras/tools"
czid_dedup_bin_path="${tools}/czid-dedup-Linux"

#main command

#check if result already exists for the current accession, if so, check gzip file integrity. If all pass, exit the script
if [ -e ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz ] && [ -e ${ddir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq.gz ] #check if dedup result files exist, PE reads
then
    if $(gzip -t ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz) && $(gzip -t ${ddir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq.gz) #check if dedup result files pass gzip integrity test, PE reads
        then
            countFileSize1=$(gzip -c ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz | wc -c | awk '{print $1}')
            countFileSize2=$(gzip -c ${ddir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq.gz | wc -c | awk '{print $1}')
            if  [ $countFileSize1 -gt 100 ] && [ $countFileSize2 -gt 100 ]
            then
                echo "Skipping ${ddir}/${accessions[$idx]} PE reads, dedup.fastq.gz found and size > 100" >> ${working_dir}/logs/dedup.skip.log
                exit 0
            fi
    fi
else #check SE
    if [ -e ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz ] #check if dedup result files exist, SE
    then
        if $(gzip -t ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz) #check if dedup result files pass gzip integrity test, SE
            then
                countFileSize1=$(gzip -c ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz | wc -c | awk '{print $1}')
                if  [ $countFileSize1 -gt 100 ]
                then
                    echo "skipping ${ddir}/${accessions[$idx]} SE reads, dedup.fastq.gz found and size > 100" >> ${working_dir}/logs/dedup.skip.log
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
if [ -e ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && [ -e ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.gz ] && [ $(gzip -l ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz | awk 'NR==2 {print $2}') -ne 0 ] && [ $(gzip -l ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.gz | awk 'NR==2 {print $2}') -ne 0 ]  #gzip file exists and is not empty
then 
    if $(gzip -t ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz) && $(gzip -t ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.gz)
    then
        echo "running dedup using ${bdir}/${accessions[$idx]} PE reads" >> ${working_dir}/logs/dedup.PE.log
        if [ ! -d "${ddir}/${accessions[$idx]}" ] 
        then
            mkdir ${ddir}/${accessions[$idx]}
        fi
        gunzip -c ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz > ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq
        gunzip -c ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.gz > ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.fastq

        #determine read length
        read_length=0
        if [ -e ${sdir}/PE/${accessions[$idx]}/Log.final.out ]
        then
            read_length_from_STAR_report=$(cat ${sdir}/PE/${accessions[$idx]}/Log.final.out | sed '7q;d'| cut -f2)
            read_length=$read_length_from_STAR_report
        else
            echo "STAR final log not found: ${sdir}/PE/${accessions[$idx]}/Log.final.out, determining read length with ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq" >> ${working_dir}/logs/dedup.warning.log
            read_length_from_fastq=$(sed '2q;d' ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq | wc -m)
            read_length=$read_length_from_fastq
        fi
        re='^[0-9]+$'
        read_length=$(echo $read_length | sed 's/[^0-9]*//g')
        if [ $read_length -eq 0 ] || ! [[ $read_length =~ $re ]] || (( $read_length < $minimum_read_length ))
        then
            echo "[Fatal error][${accessions[$idx]}] Cannot determine read length or read too short, got read length: $read_length, " >> ${working_dir}/logs/dedup.error.log
            exit 1
        fi
        length_useFor_dedup=$(bc <<< "scale=3; $(bc <<< "scale=3; $read_length / 100") * $perc_len") # calculate percent read length
        length_useFor_dedup=${length_useFor_dedup%.*}
        echo "${accessions[$idx]} read length ${read_length}, use ${length_useFor_dedup} for dedup">> ${working_dir}/logs/dedup.lengths.log

        #run dedup with -l percent read length
        ${czid_dedup_bin_path} \
        -i ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq \
        -i ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.fastq \
        -o ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq \
        -o ${ddir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq \
        -c ${ddir}/${accessions[$idx]}/cluster.csv \
        -l ${length_useFor_dedup}
        rm ${ddir}/${accessions[$idx]}/cluster.csv #  writing to cluster.csv then delete it is to preventing multiple jobs trying to write to the same file and causing unforseen problems
        rm ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq #remove unzipped BT2 results
        rm ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.fastq #remove unzipped BT2 results
        gzip ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq
        gzip ${ddir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq
        if [ $rm_BT2_res -eq 1 ] #remove bowtie2 results
        then    
            rm ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz
            rm ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.gz
        fi
    else
        echo "[Fatal error][${accessions[$idx]}] At least one of the ${bdir}/${accessions[$idx]}/*.filteredbyBT.gz files failed the integrity test" >> ${working_dir}/logs/dedup.error.log
    fi
else #SE reads
    if [ -e ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && [ $(gzip -l ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz | awk 'NR==2 {print $2}') -ne 0 ] #gzip file exists and is not empty
    then
        if $(gzip -t ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz)
        then
            echo "running dedup using ${bdir}/${accessions[$idx]} SE reads" >> ${working_dir}/logs/dedup.SE.log
            if [ ! -d "${ddir}/${accessions[$idx]}" ] 
            then
                mkdir ${ddir}/${accessions[$idx]}
            fi
            gunzip -c ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz > ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq

            #determine read length
            read_length=0
            if [ -e ${sdir}/SE/${accessions[$idx]}/Log.final.out ]
            then
                read_length_from_STAR_report=$(cat ${sdir}/SE/${accessions[$idx]}/Log.final.out | sed '7q;d'| cut -f2)
                read_length=$read_length_from_STAR_report
            else
                echo "STAR final log not found: ${sdir}/SE/${accessions[$idx]}/Log.final.out, determining read length with ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq" >> ${working_dir}/logs/dedup.warning.log
                read_length_from_fastq=$(sed '2q;d' ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq | wc -m)
                read_length=$read_length_from_fastq
            fi
            re='^[0-9]+$'
            read_length=$(echo $read_length | sed 's/[^0-9]*//g')
            if [ $read_length -eq 0 ] || ! [[ $read_length =~ $re ]] || (( $read_length < $minimum_read_length ))
            then
                echo "[Fatal error][${accessions[$idx]}] Cannot determine read length or read too short, got read length: $read_length, " >> ${working_dir}/logs/dedup.error.log
                exit 1
            fi
            length_useFor_dedup=$(bc <<< "scale=3; $(bc <<< "scale=3; $read_length / 100") * $perc_len") # calculate percent read length
            length_useFor_dedup=${length_useFor_dedup%.*}
            echo "${accessions[$idx]} read length ${read_length}, use ${length_useFor_dedup} for dedup">> ${working_dir}/logs/dedup.lengths.log

            #run dedup with -l percent read length
            ${czid_dedup_bin_path} \
            -i ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq \
            -o ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq \
            -c ${ddir}/${accessions[$idx]}/cluster.csv \
            -l ${length_useFor_dedup}
            rm ${ddir}/${accessions[$idx]}/cluster.csv # writing to cluster.csv then delete it is to preventing multiple jobs trying to write to the same file and causing unforseen problems
            rm ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.fastq #remove unzipped BT2 results
            gzip ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq
            if [ $rm_BT2_res -eq 1 ] #remove bowtie2 results
            then    rm ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz
            fi
        else
            echo "[Fatal error][${accessions[$idx]}] At least one of the ${bdir}/${accessions[$idx]}/*.filteredbyBT.gz files failed the integrity test" >> ${working_dir}/logs/dedup.error.log
        fi
    else 
        echo "[Fatal error][${accessions[$idx]}] At least one of the ${bdir}/${accessions[$idx]}/*.filteredbyBT.gz files were not found or are empty files" >> ${working_dir}/logs/dedup.error.log
    fi
fi
###########
#end dedup#
###########




