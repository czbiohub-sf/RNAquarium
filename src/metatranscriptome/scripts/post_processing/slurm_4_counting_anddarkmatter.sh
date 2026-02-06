#!/bin/bash

#SBATCH --job-name=step4_RNAquarium_metatranscriptome_darkmatteraccounting
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3200G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=16
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user eric.waltari@czbiohub.org   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc


module purge
module load anaconda/2023.03
conda activate seqtk

## generic code

## focused on dark matter calcs and extracting fasta
## /hpc/projects/balla_group/sra_experiments/RNAquarium_60k_nonhostpipeline/paired_end/chunk_001_nonzfhum.fasta


### key stats are taxonomy_hits_nonhost_inputcontigs_list
## then combined unfiltered count
## taxonomy_hits_nonhost_hitfound_list.txt
## DIFFERENCE IS DARK MATTER

# Step 1: Count sequences before filter (e.g. 89.8M contigs)
for file in chunk_[0-9][0-9][0-9].fasta; do
    grep '^>' "$file"
done > taxonomy_hits_prefilter_list.txt

wc -l taxonomy_hits_prefilter_list.txt

# Step 2: Count sequences after host filtering (e.g. 48.0M contigs)
for file in chunk_*_nonzfhum.fasta; do
    grep '^>' "$file"
done > taxonomy_hits_NTNRinput_list.txt

wc -l taxonomy_hits_NTNRinput_list.txt

## Step 3
# Clear temp file first
#> temp_nt_queries.txt

#cd nt_blast  # change to NT directory or *.txt.gz
#for file in chunk_*_nonzfhum.blast.txt.gz; do
#    zcat "$file" | awk -F'\t' '{print $1}' >> ../temp_nt_queries.txt
#done
#cd ..

#sort temp_nt_queries.txt | uniq > taxonomy_hits_blastnclustered_allhits_list.txt
#rm temp_nt_queries.txt
#wc -l taxonomy_hits_blastnclustered_allhits_list.txt


# Clean up any previous temp files
rm -f temp_nt_query_*.txt taxonomy_hits_blastnclustered_allhits_list.txt

# Loop through each NT BLAST result file
cd nt_blast
for file in chunk_*_nonzfhum.blast.txt.gz; do
    outname="../temp_nt_query_${file%.blast.txt.gz}.txt"
    zcat "$file" | awk -F'\t' '{print $1}' | sort | uniq > "$outname"
done
cd ..

# Final aggregation
cat temp_nt_query_*.txt | sort | uniq > taxonomy_hits_blastnclustered_allhits_list.txt

# Clean up
rm temp_nt_query_*.txt
wc -l taxonomy_hits_blastnclustered_allhits_list.txt


## Step 4
> temp_nr_queries.txt

cd nr_diamond  # change to NR directory was chunk_*_nonzfhum.a.diamond.txt.gz
for file in *.txt.gz; do
    zcat "$file" | awk -F'\t' '{print $1}' >> ../temp_nr_queries.txt
done
cd ..

sort temp_nr_queries.txt | uniq > taxonomy_hits_diamond_allhits_list.txt
rm temp_nr_queries.txt

wc -l taxonomy_hits_diamond_allhits_list.txt

## Step 5
cat taxonomy_hits_blastnclustered_allhits_list.txt taxonomy_hits_diamond_allhits_list.txt \
    | sort | uniq > taxonomy_hits_NTorNR_combined_list.txt

wc -l taxonomy_hits_NTorNR_combined_list.txt

## Step 6
sort taxonomy_hits_NTNRinput_list.txt > taxonomy_hits_NTNRinput_list.sorted.txt

## need to remove > at start of every line
sed 's/^>//' taxonomy_hits_NTNRinput_list.sorted.txt > taxonomy_hits_NTNRinput_list.sortedcleaned.txt

comm -23 taxonomy_hits_NTNRinput_list.sortedcleaned.txt taxonomy_hits_NTorNR_combined_list.txt \
    > taxonomy_hits_notfoundinNTorNR_list.txt

wc -l taxonomy_hits_notfoundinNTorNR_list.txt
## need count for R - taxonomy_hits_notfoundinNTorNR_count0.txt

## not quite need to then remove " taxonomy_hits_notfoundinNTorNR_list.txt"
wc -l taxonomy_hits_notfoundinNTorNR_list.txt > taxonomy_hits_notfoundinNTorNR_count.txt
sed 's/ taxonomy_hits_notfoundinNTorNR_list.txt//' taxonomy_hits_notfoundinNTorNR_count.txt > taxonomy_hits_notfoundinNTorNR_count0.txt

#comm -23 ./RNAquarium_outputs/stats/taxonomy_hits_nonhost_inputcontigs_list.txt ./RNAquarium_outputs/stats/taxonomy_hits_nonhost_hitfound_list.txt | wc -l > ./RNAquarium_outputs/stats/taxonomy_hits_notfoundinNTorNR_count0.txt



printf "Pre-filter contigs: %.1fM\n" $(($(wc -l < taxonomy_hits_prefilter_list.txt) / 1e6))
printf "Post-filter contigs: %.1fM\n" $(($(wc -l < taxonomy_hits_NTNRinput_list.txt) / 1e6))
printf "NT hits: %.1fM\n" $(($(wc -l < taxonomy_hits_blastnclustered_allhits_list.txt) / 1e6))
printf "NR hits: %.1fM\n" $(($(wc -l < taxonomy_hits_diamond_allhits_list.txt) / 1e6))
printf "NT or NR hits (union): %.1fM\n" $(($(wc -l < taxonomy_hits_NTorNR_combined_list.txt) / 1e6))
printf "Dark matter contigs: %.1fM\n" $(($(wc -l < taxonomy_hits_notfoundinNTorNR_list.txt) / 1e6))

pre=$(wc -l < taxonomy_hits_prefilter_list.txt)
post=$(wc -l < taxonomy_hits_NTNRinput_list.txt)
nt=$(wc -l < taxonomy_hits_blastnclustered_allhits_list.txt)
nr=$(wc -l < taxonomy_hits_diamond_allhits_list.txt)
union=$(wc -l < taxonomy_hits_NTorNR_combined_list.txt)
dark=$(wc -l < taxonomy_hits_notfoundinNTorNR_list.txt)

# Write values to file for Python to read
echo "$pre $post $nt $nr $union $dark" > pipeline_counts.txt


## this loop will find all of the fasta files in the work directory, then run seqtk
## TO GET THESE DARK MATTER FASTAS, WE CAN USE NONZFHUM FASTAS (in paired_end folder)

#find ./paired_end -type f -name "*nonzfhum.fasta"
#find ./paired_end -type f -name "*hum.fasta" | while read -r fasta_file; do
#    seqtk subseq "$fasta_file" ./RNAquarium_outputs/stats/taxonomy_hits_notfoundinNTorNR_list.txt >> #./RNAquarium_outputs/taxonomy_hits_notfoundinNTorNR_darkmatter.fasta
#done


for file in chunk_*_nonzfhum.fasta; do
    seqtk subseq "$file" taxonomy_hits_notfoundinNTorNR_list.txt >> taxonomy_hits_notfoundinNTorNR_darkmatter.fasta
done

### final ideas - to get the filtered list per killifish code, then create a simple workflow figure incorporating all of these numbers using pikchr! would also need count of contigs input into first blast...

## can add plotting steps via R commands invoked last


sleep 5
module purge
module load R/4.4
#R
#source("nr_diamond_processing.r")
Rscript 4_alluvial_withdarkmatter.r

sleep 5
cd RNAquarium_outputs
mkdir -p stats
cd ..

mv taxonomy_hits_prefilter_list.txt RNAquarium_outputs/stats/
mv taxonomy_hits_NTNRinput_list.txt RNAquarium_outputs/stats/
mv taxonomy_hits_blastnclustered_allhits_list.txt RNAquarium_outputs/stats/
mv taxonomy_hits_diamond_allhits_list.txt RNAquarium_outputs/stats/
mv taxonomy_hits_NTorNR_combined_list.txt RNAquarium_outputs/stats/
mv taxonomy_hits_notfoundinNTorNR_count.txt RNAquarium_outputs/stats/
mv taxonomy_hits_notfoundinNTorNR_count0.txt RNAquarium_outputs/stats/

mv taxonomy_hits_NTNRinput_list.sorted.txt RNAquarium_outputs/stats/
mv taxonomy_hits_NTNRinput_list.sortedcleaned.txt RNAquarium_outputs/stats/
mv taxonomy_hits_notfoundinNTorNR_list.txt RNAquarium_outputs/stats/
mv pipeline_counts.txt RNAquarium_outputs/stats/

#mv taxonomy_hits_notfoundinNTorNR_list.txt stats/
mv taxonomy_hits_notfoundinNTorNR_darkmatter.fasta RNAquarium_outputs/
