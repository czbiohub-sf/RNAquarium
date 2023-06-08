#!/bin/bash

#SBATCH --job-name=sra_gsnap
#SBATCH --time=14-00:00:00
#SBATCH --array=1-526%500
#SBATCH --nodes=1
#SBATCH --partition preempted
#SBATCH --ntasks=1
#SBATCH --mem=96G
#SBATCH --cpus-per-task=2
#SBATCH -e slurm.out.gsnap/slurm-%A_%a.err
#SBATCH -o slurm.out.gsnap/slurm-%A_%a.out

declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

module load anaconda
conda activate zf_pipeline

#parameters
max_mismatch=0.3
rm_dedup_res=0 #whether to remove dedup results, 0=don't remove, 1=remove

#setting directories
#working_dir="/hpc/projects/balla_group/sra_experiments/all_zebrafish_RNAseq/unmapped_dev"
working_dir="/hpc/scratch/group.theory/jparas/zf_pipeline"
ddir=${working_dir}/Dedup_out
gdir=${working_dir}/Gsnap_out

# declare arrays
#readarray -t accessions < <(cat /hpc/projects/balla_group/sra_experiments/all_zebrafish_RNAseq/SRA_accession_list.1.27.23.txt)
readarray -t accessions < <(cat "${working_dir}/data/gsnap.missing.txt")

tools="/hpc/projects/theory/sharing/internship/jacob.paras/tools"
bbmap_dir="${tools}/bbmap"
gsnap_out_bin_path="${tools}/gmap-2021-12-17/bin/gsnap" #For small genomes of less than 2^32 (4 billion) bp, please run gsnap
gdbdir="${tools}/gmap-2021-12-17/db"

gdbidx="Danio_rerio.GRCz11.dna_sm.primary_assembly"

#check if gsnap result already exists for the current accession, if so, check gzip file integrity. If all pass, exit the script
# if [ -e ${gdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz ] && [ -e ${gdir}/${accessions[$idx]}/gsnap.stats.txt ] && [ -e ${gdir}/${accessions[$idx]}/gsnap.mapped.names.txt ] 
# then
#     if $(gzip -t ${gdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz) 
#     then
#         countFileSize=$(gzip -c ${gdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz | wc -c | awk '{print $1}')
#         countFileSize2=$(wc -c ${gdir}/${accessions[$idx]}/gsnap.stats.txt | awk '{print $1}')
#         countFileSize3=$(wc -c ${gdir}/${accessions[$idx]}/gsnap.mapped.names.txt | awk '{print $1}')
#         if [ $countFileSize -gt 100 ] && [ $countFileSize2 -gt 4 ] && [ $countFileSize3 -gt 4 ]
#         then
#             echo "gsnap result already exists for ${accessions[$idx]} and gzip file integrity and sizes are good, skipping gsnap" >> ${working_dir}/logs/gsnap.skip.log
#             exit 0
#         fi
#     fi
# fi


# PE reads
if [ -e ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz ] && [ -e ${ddir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq.gz ]
then 
    if $(gzip -t ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz) && $(gzip -t ${ddir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq.gz)
    then
        echo "running gsnap using ${ddir}/${accessions[$idx]} PE reads" >> ${working_dir}/logs/gsnap.process.log
        if [ ! -d "${gdir}/${accessions[$idx]}" ] 
        then
            mkdir ${gdir}/${accessions[$idx]}
        fi
        #gsnap command
        ${gsnap_out_bin_path} \
            -A sam \
            -N 1 `#find novel splice sites, this optin will turn on the RNAseq mode without the need to supply known splice sites` \
            --batch=2 `# tried mode 4, and several runs returned OOM errors` \
            --use-shared-memory=0 `#private memory for each instance` \
            --npaths=1 `#Maximum number of paths to print (default 100).` \
            --ordered `#Print output in same order as input (relevant only if there is more than one worker thread)` \
            -t 2 `#thread` \
            --max-mismatches=${max_mismatch} `#czID uses 40bp ( defult is 30% of read length)` \
            -D ${gdbdir}  `#Genome directory.  Default (as specified by --with-gmapdb to the configure program) is /hpc/projects/balla_group/sra_experiments/tools/gmap-2021-12-17/db` \
            -d ${gdbidx} \
            -o ${gdir}/${accessions[$idx]}/gsnap_out.sam \
            --gunzip `#Uncompress gzipped input files` \
            ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz ${ddir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq.gz  `# Two input FASTAs means paired reads.` \
        #process Sam file
            if [ -e ${gdir}/${accessions[$idx]}/gsnap_out.sam ]
            then #get stats
                samtools view ${gdir}/${accessions[$idx]}/gsnap_out.sam | cut -f2 | sort | uniq -c > ${gdir}/${accessions[$idx]}/gsnap.stats.txt
                # 0: mapped forward (unpaired)
                # 16: mapped reverse (unpaired)
                # 4: unmapped (unpaired)
                #get name of mapped reads
                samtools view ${gdir}/${accessions[$idx]}/gsnap_out.sam | awk '{if( and($2,4)==0 && and($2,8)==0) {print $1}}' > ${gdir}/${accessions[$idx]}/gsnap.mapped.names.txt # print out read name if both read and mate are mapped ( that is SAM bitwise flag does not contain 4 or 8 )

                #filter dedup out using names of mapped reads
                gunzip -c ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz | ${bbmap_dir}/filterbyname.sh in=stdin names=${gdir}/${accessions[$idx]}/gsnap.mapped.names.txt out=${gdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq include=f
                gunzip -c ${ddir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq.gz | ${bbmap_dir}/filterbyname.sh in=stdin names=${gdir}/${accessions[$idx]}/gsnap.mapped.names.txt out=${gdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.gsnapFiltered.fastq include=f
                #compress filtered reads
                gzip ${gdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq
                gzip ${gdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.gsnapFiltered.fastq

                #remove dedup result out
                if [ $rm_dedup_res -eq 1 ]
                then
                    if [ -e ${gdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz ] && $(gzip -t ${gdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz)
                    then
                        rm ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz
                    fi
                    if [ -e ${gdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.gsnapFiltered.fastq.gz ] && $(gzip -t ${gdir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.gsnapFiltered.fastq.gz)
                    then
                        rm ${ddir}/${accessions[$idx]}/Unmapped.out.mate2.filteredbyBT.dedup.fastq.gz
                    fi
                fi

                #remove gsnap bam file
                rm ${gdir}/${accessions[$idx]}/gsnap_out.sam
            fi
    else
        echo "At least one of the ${ddir}/${accessions[$idx]}/*.filteredbyBT.dedup.fastq.gz files failed the integrity test" >> ${working_dir}/logs/gsnap.error.log
    fi
else #SE reads
    if [ -e ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz ]
    then
        if $(gzip -t ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz)
        then
            echo "running gsnap using ${ddir}/${accessions[$idx]} SE reads" >> ${working_dir}/logs/gsnap.process.log
            if [ ! -d "${gdir}/${accessions[$idx]}" ] 
            then
                mkdir ${gdir}/${accessions[$idx]}
            fi
        #gsnap command
        ${gsnap_out_bin_path} \
            -A sam \
            -N 1 `#find novel splice sites, this optin will turn on the RNAseq mode without the need to supply known splice sites` \
            --batch=2 `# tried mode 4, and several runs returned OOM errors` \
            --use-shared-memory=0 `#private memory for each instance` \
            --npaths=1 `#Maximum number of paths to print (default 100).` \
            --ordered `#Print output in same order as input (relevant only if there is more than one worker thread)` \
            -t 2 `#thread` \
            --max-mismatches=${max_mismatch} `#czID uses 40bp ( defult is 30% of read length)` \
            -D ${gdbdir}  `#Genome directory.  Default (as specified by --with-gmapdb to the configure program) is /hpc/projects/balla_group/sra_experiments/tools/gmap-2021-12-17/db` \
            -d ${gdbidx} \
            -o ${gdir}/${accessions[$idx]}/gsnap_out.sam \
            --gunzip `#Uncompress gzipped input files` \
            ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz `# Two input FASTAs means paired reads.` \
        #process Sam file
            if [ -e ${gdir}/${accessions[$idx]}/gsnap_out.sam ]
            then #get stats
                samtools view ${gdir}/${accessions[$idx]}/gsnap_out.sam | cut -f2 | sort | uniq -c > ${gdir}/${accessions[$idx]}/gsnap.stats.txt
                # 0: mapped forward (unpaired)
                # 16: mapped reverse (unpaired)
                # 4: unmapped (unpaired)
                #get name of mapped reads
                samtools view ${gdir}/${accessions[$idx]}/gsnap_out.sam | awk '{if( and($2,4)==0 && and($2,8)==0) {print $1}}' > ${gdir}/${accessions[$idx]}/gsnap.mapped.names.txt # print out read name if both read and mate are mapped ( that is SAM bitwise flag does not contain 4 or 8 )
                #filter dedup out using names of mapped reads
                gunzip -c ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz | ${bbmap_dir}/filterbyname.sh in=stdin names=${gdir}/${accessions[$idx]}/gsnap.mapped.names.txt out=${gdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq include=f
                #compress filtered reads
                gzip ${gdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq

                #remove dedup result out
                if [ $rm_dedup_res -eq 1 ]
                then
                    if [ -e ${gdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz ] && $(gzip -t ${gdir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz)
                    then
                        rm ${ddir}/${accessions[$idx]}/Unmapped.out.mate1.filteredbyBT.dedup.fastq.gz
                    fi
                fi

                #remove gsnap bam file
                rm ${gdir}/${accessions[$idx]}/gsnap_out.sam
            fi
        else
            echo "At least one of the ${ddir}/${accessions[$idx]}/*.filteredbyBT.dedup.fastq.gz files failed the integrity test" >> ${working_dir}/logs/gsnap.error.log
        fi
    else 
        echo "At least one of the ${ddir}/${accessions[$idx]}/*.filteredbyBT.dedup.fastq.gz files were not found" >> ${working_dir}/logs/gsnap.error.log
    fi
fi


#gsnap --batch parameter:
#   -B, --batch=INT                Batch mode (default = 2)
#                                  Mode     Hash offsets  Hash positions  Genome          Local hash offsets  Local hash positions
#                                    0      allocate      mmap            mmap            allocate            mmap
#                                    1      allocate      mmap & preload  mmap            allocate            mmap & preload
#                                    2      allocate      mmap & preload  mmap & preload  allocate            mmap & preload
#                                    3      allocate      allocate        mmap & preload  allocate            allocate
#                       (default)    4      allocate      allocate        allocate        allocate            allocate




