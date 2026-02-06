#!/bin/bash

#SBATCH --job-name=step1b_RNAquarium_metatranscriptome_NTNRcombining_andplotting_afterbbdukfilter
#SBATCH --time=22:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3900G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=32
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user eric.waltari@czbiohub.org   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc

## just the bash commands I would use

module purge
module load R/4.4
#mkdir RNAquarium_outputs
#R
#source("nr_diamond_processing.r")

echo "Starting R script at $(date)"
Rscript 1_combinecontigs_alluvial_treemap_bbdukfilter.r

cd RNAquarium_outputs

# Validate input files taxonomy_hits_viruses0_fullcols_mostrecent.tsv
if [ ! -f taxonomy_hits_viruses0_fullcols_mostrecent.tsv ]; then
    echo "Error: TSV file from R not found: taxonomy_hits_viruses0_fullcols_mostrecent.tsv"
    exit 1
fi


#gzip -k taxonomy_hits_viruses0_fullcols_mostrecent.tsv
gzip -f taxonomy_hits_viruses0_fullcols_mostrecent.tsv

cd ..
chmod +x taxname_to_taxid_viruses0_twoparts.sh
chmod +x taxname_to_taxid_nonhost_twoparts.sh



echo "Running taxname to taxid script for viruses at $(date)"
./taxname_to_taxid_viruses0_twoparts.sh

sleep 1
echo "Running taxname to taxid script for all non-host at $(date)"
./taxname_to_taxid_nonhost_twoparts.sh


sleep 3
cd RNAquarium_outputs

## renames and unzips
## rename original to _pretaxid taxonomy_hits_nonhost_mostrecent.tsv.gz
mv taxonomy_hits_viruses0_fullcols_mostrecent.tsv.gz taxonomy_hits_viruses0_fullcols_mostrecent_pretaxid.tsv.gz
#mv taxonomy_hits_viruses0_fullcols_mostrecent.tsv taxonomy_hits_viruses0_fullcols_mostrecent_pretaxid.tsv
mv taxonomy_hits_nonhost_mostrecent.tsv.gz taxonomy_hits_nonhost_mostrecent_pretaxid.tsv.gz

# rename taxonomy_hits_viruses0_fullcols_mostrecent_withtaxids.tsv.gz to taxonomy_hits_viruses0_fullcols_mostrecent_withtaxids_part1.tsv.gz

mv taxonomy_hits_viruses0_fullcols_mostrecent_withtaxids.tsv.gz taxonomy_hits_viruses0_fullcols_mostrecent_withtaxids_part1.tsv.gz
mv taxonomy_hits_viruses0_fullcols_mostrecent_withtaxids_final.tsv.gz taxonomy_hits_viruses0_fullcols_mostrecent.tsv.gz

mv taxonomy_hits_nonhost_mostrecent_withtaxids.tsv.gz taxonomy_hits_nonhost_mostrecent_withtaxids_part1.tsv.gz
mv taxonomy_hits_nonhost_mostrecent_withtaxids_final.tsv.gz taxonomy_hits_nonhost_mostrecent.tsv.gz


# rename to og taxonomy_hits_viruses0_fullcols_mostrecent_withtaxids_final.tsv.gz
## then unzip 
gunzip -k -f taxonomy_hits_viruses0_fullcols_mostrecent.tsv.gz
sleep 3
## REMEMBER IN LAST R SCRIPT TO START BY LOADING UPDATED FILES

cd ..
## then last R script
echo "Starting final R script at $(date)"
Rscript 1_combinecontigs_alluvial_treemap_withtaxid_bbdukfilter.r

## now new commands related to taxname_to_taxid script!
## need to first gzip virus0 file
## also at end of scripts rename both OG tsv files and then rename newones to original names

## finally run last bit of r script
## also repeat all of these steps for slurm_step1 as well...

echo "All done at $(date)"