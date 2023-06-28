#!/bin/bash

#SBATCH --job-name=sra_bowtie
#SBATCH --time=14-00:00:00
#SBATCH --array=1-380%400
#SBATCH --partition preempted
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=2
#SBATCH -e slurm.out.Bowtie2.missing/slurm-%A_%a.err
#SBATCH -o slurm.out.Bowtie2.missing/slurm-%A_%a.out

declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

module load anaconda
conda activate zf_pipeline

#use /tmp
use_tmp=0

#remove STAR_out after Bowtie2
remove_STAR_OUT=0

#setting directories
working_dir="/hpc/scratch/group.theory/jparas/zf_pipeline"
fdir=${working_dir}/fastq
sdir=${working_dir}/STAR_out
bdir=${working_dir}/Bowtie2_out

# declare arrays
readarray -t accessions < <(cat "${working_dir}/data/Bowtie2.missing.txt")

tools="/hpc/projects/theory_ds/internship/jacob.paras/tools"
bbmap_dir="${tools}/bbmap"
samtools="${tools}/samtools-1.6/samtools"

#setting up index
ZF_idx="${tools}/bowtie2_index/Danio_rerio.GRCz11.dna_sm.primary_assembly"


####################################################
# check if Bowtie2 exists, skip current job if so #
####################################################
# # PE
# if [ -e ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && [ -e ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.gz ] && [ -e ${bdir}/${accessions[$idx]}/bowtie2.stats.txt ]
# then 
#     if $(gzip -t ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz) && $(gzip -t ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.gz)
#     then
#         countFileSize=$(gzip -c ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz | wc -c | awk '{print $1}')
#         if [ $countFileSize -gt 100 ]
#         then
#             echo "${accessions[$idx]} Bowtie2_out PE result exists and size > 100" >> ${working_dir}/logs/Bowtie2.skip.log
#             exit 0
#         fi
#     fi
# fi
# # SE
# if [ -e ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && [ -e ${bdir}/${accessions[$idx]}/bowtie2.stats.txt ]
# then 
#     if $(gzip -t ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz) 
#     then
#         countFileSize=$(gzip -c ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz | wc -c | awk '{print $1}')
#         if [ $countFileSize -gt 100 ]
#         then
#             echo "${accessions[$idx]} Bowtie2_out SE result exists and size > 100" >> ${working_dir}/logs/Bowtie2.skip.log
#             exit 0
#         fi
#     fi
# fi

###########################
# start executing Bowtie2 #
###########################

if [ $use_tmp -eq 1 ]
then
    bdir=/tmp/Bowtie2_out

    if [ ! -d $bdir ]
    then
        mkdir $bdir
    fi

fi


if [ -d ${sdir}/PE/${accessions[$idx]} ]
then  #paired end
    echo "running bowtie2 with ${accessions[$idx]} PE STAR_out" >> ${working_dir}/logs/Bowtie2.PE.log
    #make bowtie output dir
    rm -rf ${bdir}/${accessions[$idx]}
    mkdir ${bdir}/${accessions[$idx]}
    #map mate 1 (in unpaired mode)
    if [ -e ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate1.gz ]
    then 
        echo "running bowtie2 with ${accessions[$idx]} PE STAR_out Unmapped.out.mate1.gz"
        bowtie2 --quiet --very-sensitive-local \
        --rg-id na --rg LB:na --rg SM:na --rg PL:na --rg PU:na \
        -x ${ZF_idx} -U ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate1.gz | \
        $samtools view -b - >> ${bdir}/${accessions[$idx]}/bowtie2.bam
    fi
    #map mate 2 (in unpaired mode)
    if [ -e ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate2.gz ]
    then 
        echo "running bowtie2 with ${accessions[$idx]} PE STAR_out Unmapped.out.mate2.gz"
        bowtie2 --quiet --very-sensitive-local \
        --rg-id na --rg LB:na --rg SM:na --rg PL:na --rg PU:na \
        -x ${ZF_idx} -U ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate2.gz | \
        $samtools view -b - >> ${bdir}/${accessions[$idx]}/bowtie2.bam
    fi
    #process bam file
    if [ -e ${bdir}/${accessions[$idx]}/bowtie2.bam ]
    then #get stats
        echo "processing bowtie2 bam file"
        $samtools view ${bdir}/${accessions[$idx]}/bowtie2.bam | cut -f2 | sort | uniq -c > ${bdir}/${accessions[$idx]}/bowtie2.stats.txt
        # retrain unmapped reads if bitwise flag contains 4 or 8
        # 4: unmapped 
        # 8: mate unmapped
        $samtools view ${bdir}/${accessions[$idx]}/bowtie2.bam | awk '{if( and($2,4)==0 && and($2,8)==0) {print $1}}' > ${bdir}/${accessions[$idx]}/bowtie2.mapped.names.txt # print out read name if both read and mate are mapped ( that is SAM bitwise flag does not contain 4 or 8 )
        #filter STAR out using names of mapped reads
        gunzip -c ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate1.gz | ${bbmap_dir}/filterbyname.sh in=stdin names=${bdir}/${accessions[$idx]}/bowtie2.mapped.names.txt out=${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT include=f
        gunzip -c ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate2.gz | ${bbmap_dir}/filterbyname.sh in=stdin names=${bdir}/${accessions[$idx]}/bowtie2.mapped.names.txt out=${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT include=f
        #compress filtered reads
        gzip ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT
        gzip ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT

        #remove STAR out
        if [ $remove_STAR_OUT -eq 1 ]
        then
            if [ -e ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && $(gzip -t ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz)
            then
                rm ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate1.gz
            fi
            if [ -e ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.gz ] && $(gzip -t ${bdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.gz)
            then
                rm ${sdir}/PE/${accessions[$idx]}/Unmapped.out.mate2.gz
            fi
        fi
    #remove bowtie bam
    rm ${bdir}/${accessions[$idx]}/bowtie2.bam

    fi
elif [ -d ${sdir}/SE/${accessions[$idx]} ]
    then #single end
        echo "running bowtie2 with ${accessions[$idx]} SE STAR_out found" >> ${working_dir}/logs/Bowtie2.SE.log
        #make bowtie output dir
        rm -rf ${bdir}/${accessions[$idx]}
        mkdir ${bdir}/${accessions[$idx]}

        #map in unpaired mode
        if [ -e ${sdir}/SE/${accessions[$idx]}/Unmapped.out.mate1.gz ]
        then 
            echo "running bowtie2 with ${accessions[$idx]} SE STAR_out Unmapped.out.mate1.gz"
            bowtie2 --quiet --very-sensitive-local \
            --rg-id na --rg LB:na --rg SM:na --rg PL:na --rg PU:na \
            -x ${ZF_idx} -U ${sdir}/SE/${accessions[$idx]}/Unmapped.out.mate1.gz | \
            $samtools view -b - >> ${bdir}/${accessions[$idx]}/bowtie2.bam
        fi
        
        #process bam file
        if [ -e ${bdir}/${accessions[$idx]}/bowtie2.bam ]
        then #get stats
            echo "processing bowtie2 bam file"
            $samtools view ${bdir}/${accessions[$idx]}/bowtie2.bam | cut -f2 | sort | uniq -c > ${bdir}/${accessions[$idx]}/bowtie2.stats.txt
            # retain unmapped reads if bitwise flag contains 4 or 8
            # 4: unmapped 
            # 8: mate unmapped
            $samtools view ${bdir}/${accessions[$idx]}/bowtie2.bam | awk '{if( and($2,4)==0 && and($2,8)==0) {print $1}}' > ${bdir}/${accessions[$idx]}/bowtie2.mapped.names.txt # print out read name if both read and mate are mapped ( that is SAM bitwise flag does not contain 4 or 8 )
            #filter STAR out using names of mapped reads
            gunzip -c ${sdir}/SE/${accessions[$idx]}/Unmapped.out.mate1.gz | ${bbmap_dir}/filterbyname.sh in=stdin names=${bdir}/${accessions[$idx]}/bowtie2.mapped.names.txt out=${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT include=f
            
            #compress filtered reads
            gzip ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT

            #remove STAR out
            if [ $remove_STAR_OUT -eq 1 ]
            then
                if [ -e ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz ] && $(gzip -t ${bdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.gz)
                then
                    rm ${sdir}/SE/${accessions[$idx]}/Unmapped.out.mate1.gz
                fi
            fi
        fi
        #remove bowtie bam
        rm ${bdir}/${accessions[$idx]}/bowtie2.bam

else
    echo "${accessions[$idx]} was not found in STAR_out" >> ${working_dir}/logs/Bowtie2.fail.log
    echo "${accessions[$idx]} was not found in STAR_out"
fi


# if using /tmp, copy result files out of /tmp
if [ ${use_tmp} -eq 1 ]
then
    bdir2=${working_dir}/Bowtie2_out

    if [ -e "${bdir}/${accessions[$idx]}" ]
    then
        echo "copying ${accessions[$idx]} from /tmp to ${bdir2}"
        rm -rf ${bdir2}/${accessions[$idx]}
        mv ${bdir}/${accessions[$idx]} ${bdir2}/${accessions[$idx]}
    fi

fi