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

## adding checks to see if parts 1-4 already exist!
# Step 1–4 checkpoint file names
pre_file="taxonomy_hits_prefilter_list.txt"
post_file="taxonomy_hits_NTNRinput_list.txt"
nt_file="taxonomy_hits_blastnclustered_allhits_list.txt"
nr_file="taxonomy_hits_diamond_allhits_list.txt"

## October update to be more thorough
# Check if all step 1–4 outputs already exist
# Enable safe globbing (no literal patterns if no match)
shopt -s nullglob

if [[ -f "$pre_file" && -f "$post_file" && -f "$nt_file" && -f "$nr_file" ]]; then
    echo "✔ Step 1–4 output files already exist. Skipping to step 5..."
else
    echo "⏳ Step 1–4 outputs not found. Running initial steps..."

    ##########################
    # Step 1: Count sequences before host filtering
    echo "📍 Step 1: Counting pre-filter contigs..."

    fasta_files=(chunk_[0-9][0-9][0-9].fasta)
    if [ ${#fasta_files[@]} -eq 0 ]; then
        echo "⚠️  WARNING: No pre-filter FASTA files found matching chunk_[0-9][0-9][0-9].fasta"
    else
        for file in "${fasta_files[@]}"; do
            grep '^>' "$file"
        done > taxonomy_hits_prefilter_list.txt
        wc -l taxonomy_hits_prefilter_list.txt
    fi

    ##########################
    # Step 2: Count sequences after host filtering
    echo "📍 Step 2: Counting post-filter contigs..."

    filtered_files=(chunk_*_nonzfhum.fasta)
    if [ ${#filtered_files[@]} -eq 0 ]; then
        echo "⚠️  WARNING: No filtered FASTA files found matching chunk_*_nonzfhum.fasta"
    else
        for file in "${filtered_files[@]}"; do
            grep '^>' "$file"
        done > taxonomy_hits_NTNRinput_list.txt
        wc -l taxonomy_hits_NTNRinput_list.txt
    fi

    ##########################
    # Step 3: Process NT BLAST hits
    echo "📍 Step 3: Processing NT BLAST hits..."
    rm -f temp_nt_query_*.txt
    cd nt_blast || { echo "❌ ERROR: Could not enter nt_blast directory"; exit 1; }

    blast_files=(chunk_*_nonzfhum.blast.txt.gz)
    if [ ${#blast_files[@]} -eq 0 ]; then
        echo "⚠️  WARNING: No BLAST files found matching chunk_*_nonzfhum.blast.txt.gz"
    else
        for file in "${blast_files[@]}"; do
            outname="../temp_nt_query_${file%.blast.txt.gz}.txt"
            zcat "$file" | awk -F'\t' '{print $1}' | sort | uniq > "$outname"
        done
    fi
    cd ..

    if ls temp_nt_query_*.txt 1> /dev/null 2>&1; then
        cat temp_nt_query_*.txt | sort | uniq > taxonomy_hits_blastnclustered_allhits_list.txt
        rm temp_nt_query_*.txt
        wc -l taxonomy_hits_blastnclustered_allhits_list.txt
    else
        echo "⚠️  WARNING: No temp NT query files were created. Skipping BLAST aggregation."
    fi

    ##########################
    # Step 4: Process NR DIAMOND hits
    echo "📍 Step 4: Processing NR DIAMOND hits..."
    > temp_nr_queries.txt
    cd nr_diamond || { echo "❌ ERROR: Could not enter nr_diamond directory"; exit 1; }

    diamond_files=(*.txt.gz)
    if [ ${#diamond_files[@]} -eq 0 ]; then
        echo "⚠️  WARNING: No DIAMOND result files (*.txt.gz) found in nr_diamond/"
    else
        for file in "${diamond_files[@]}"; do
            zcat "$file" | awk -F'\t' '{print $1}' >> ../temp_nr_queries.txt
        done
    fi
    cd ..

    if [ -s temp_nr_queries.txt ]; then
        sort temp_nr_queries.txt | uniq > taxonomy_hits_diamond_allhits_list.txt
        rm temp_nr_queries.txt
        wc -l taxonomy_hits_diamond_allhits_list.txt
    else
        echo "⚠️  WARNING: temp_nr_queries.txt is empty. Skipping DIAMOND aggregation."
    fi

    echo "✅ Steps 1–4 completed (with warnings if shown above)."
fi

## Step 5
cat taxonomy_hits_blastnclustered_allhits_list.txt taxonomy_hits_diamond_allhits_list.txt \
    | sort | uniq > taxonomy_hits_NTorNR_combined_list.txt

wc -l taxonomy_hits_NTorNR_combined_list.txt


## Step 5b?
## need to take into account chunks_diamond_hits_fishorprimatescombined.txt and chunks_blastnclustered_hits_fishorprimatescombined.txt
## actually no - these are accounted for by using raw blast & diamond outputs for taxonomy_hits_NTorNR_combined_list
# cat chunks_blastnclustered_hits_fishorprimatescombined.txt chunks_diamond_hits_fishorprimatescombined.txt \
#     | sort | uniq > taxonomy_hits_NTorNR_combined_fishorprimateslist0.txt
# ## remove any from fishorprimateslist that are also in taxonomy_hits_NTorNR_combined_list
# ## check in ntblast chunk001 for PRJNA630877_S_NODE_112861
# grep "PRJNA630877_S_NODE_112861" chunk_001_blastnclustered_hits_nofishnohuman_2025-10-04.tab ## nothing found
# ## checking raw blast output???
# zcat chunk_001_nonzfhum.blast.txt.gz | grep "PRJNA630877_S_NODE_112861"
# ## it is here, it's another fish non-Danio
# ## BUT IT IS STILL HERE
# grep "PRJNA630877_S_NODE_112861" taxonomy_hits_NTorNR_combined_list.txt
# ## how about here in full nt taxonomy_hits_blastnclustered_allhits_list.txt
# grep "PRJNA630877_S_NODE_112861" taxonomy_hits_blastnclustered_allhits_list.txt
# comm -23 taxonomy_hits_NTorNR_combined_fishorprimateslist0.txt taxonomy_hits_NTorNR_combined_list.txt \
#     > taxonomy_hits_NTorNR_combined_fishorprimateslist.txt

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


## also note taxonomy_hits_nonhost_sequences_to_exclude.txt
cp RNAquarium_outputs/taxonomy_hits_nonhost_sequences_to_exclude.txt taxonomy_hits_nonhost_sequences_to_exclude.txt
## don't forget to sort!
sort taxonomy_hits_nonhost_sequences_to_exclude.txt > taxonomy_hits_nonhost_sequences_to_exclude_sorted.txt
## remove the exclude list to get the count post-BBDuk filter
comm -23 taxonomy_hits_NTorNR_combined_list.txt taxonomy_hits_nonhost_sequences_to_exclude_sorted.txt > taxonomy_hits_NTorNR_combined_list_postfilter.txt


printf "Pre-filter contigs: %.1fM\n" $(($(wc -l < taxonomy_hits_prefilter_list.txt) / 1e6))
printf "Post-filter contigs: %.1fM\n" $(($(wc -l < taxonomy_hits_NTNRinput_list.txt) / 1e6))
printf "NT hits: %.1fM\n" $(($(wc -l < taxonomy_hits_blastnclustered_allhits_list.txt) / 1e6))
printf "NR hits: %.1fM\n" $(($(wc -l < taxonomy_hits_diamond_allhits_list.txt) / 1e6))
printf "NT or NR hits (union): %.1fM\n" $(($(wc -l < taxonomy_hits_NTorNR_combined_list.txt) / 1e6))
printf "Contigs removed after BBduk filter & min150bp cutoff: %.1fM\n" $(($(wc -l < taxonomy_hits_nonhost_sequences_to_exclude.txt) / 1e6))
printf "Contigs remaining after BBduk filter & min150bp cutoff: %.1fM\n" $(($(wc -l < taxonomy_hits_NTorNR_combined_list_postfilter.txt) / 1e6))
printf "Dark matter contigs: %.1fM\n" $(($(wc -l < taxonomy_hits_notfoundinNTorNR_list.txt) / 1e6))

pre=$(wc -l < taxonomy_hits_prefilter_list.txt)
post=$(wc -l < taxonomy_hits_NTNRinput_list.txt)
nt=$(wc -l < taxonomy_hits_blastnclustered_allhits_list.txt)
nr=$(wc -l < taxonomy_hits_diamond_allhits_list.txt)
union=$(wc -l < taxonomy_hits_NTorNR_combined_list.txt)
bbdukfiltered=$(wc -l < taxonomy_hits_NTorNR_combined_list_postfilter.txt)
dark=$(wc -l < taxonomy_hits_notfoundinNTorNR_list.txt)

# Write values to file for Python to read
echo "$pre $post $nt $nr $union $bbdukfiltered $dark" > pipeline_counts.txt


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
Rscript 4b_alluvial_withdarkmatter.r

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

mv taxonomy_hits_NTorNR_combined_list_postfilter.txt RNAquarium_outputs/stats/
#mv taxonomy_hits_nonhost_sequences_to_exclude_sorted.txt RNAquarium_outputs/stats/

#mv taxonomy_hits_notfoundinNTorNR_list.txt stats/
mv taxonomy_hits_notfoundinNTorNR_darkmatter.fasta RNAquarium_outputs/
## only in 4b
mv taxonomy_hits_nonhost_sequences_to_exclude_sorted.txt RNAquarium_outputs/stats/