#!/bin/bash

#SBATCH --job-name=stepv10_RNAquarium_metatranscriptome_viralclusteringminimap
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=130G
#SBATCH --partition cpu         # Partition to submit
#SBATCH --gpus 0            # Reserve 0 GPUs for usage
#SBATCH --cpus-per-task=8
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user eric.waltari@czbiohub.org   # this is the email you wish to be notified at
#SBATCH --mail-type END,FAIL   # ALL will alert you of job beginning, completion, failure etc

## just the bash commands I would use


## then finally minimap2 will be yet separate...

module purge
module load minimap2/2.26


### then loops to pull first sequence & run minimap
## start in usual folder

#mv ./RNAquarium_outputs/virus_outputs/clusters_forminimap/taxonomy_hits_viruses_clustered_singletons* ./RNAquarium_outputs/virus_outputs


for file in ./RNAquarium_outputs/virus_outputs/clusters_forminimap/*.fasta; do
    base_name=$(basename "$file" .fasta)
    output_dir="./RNAquarium_outputs/virus_outputs/clusters_forminimap"  # Using the same base path dynamically
    mkdir -p "$output_dir"  # Ensure the output directory exists
    head -n 2 "$file" > "$output_dir/${base_name}.fa"
done


## now minimap2 command for each

for file in ./RNAquarium_outputs/virus_outputs/clusters_forminimap/*.fasta; do
    base_name=$(basename "$file" .fasta)
    output_dir="./RNAquarium_outputs/virus_outputs/clusters_forminimap"  # Using the same base path dynamically
    mkdir -p "$output_dir"  # Ensure the output directory exists
    minimap2 -ax map-ont "$output_dir/${base_name}.fa" "$file" > "$output_dir/${base_name}.sam"
done



#####################################
## need to repeat for phage subfolder

for file in ./RNAquarium_outputs/virus_outputs/clusters_forminimap_phage/*.fasta; do
    base_name=$(basename "$file" .fasta)
    output_dir="./RNAquarium_outputs/virus_outputs/clusters_forminimap_phage"  # Using the same base path dynamically
    mkdir -p "$output_dir"  # Ensure the output directory exists
    head -n 2 "$file" > "$output_dir/${base_name}.fa"
done


## now minimap2 command for each

for file in ./RNAquarium_outputs/virus_outputs/clusters_forminimap_phage/*.fasta; do
    base_name=$(basename "$file" .fasta)
    output_dir="./RNAquarium_outputs/virus_outputs/clusters_forminimap_phage"  # Using the same base path dynamically
    mkdir -p "$output_dir"  # Ensure the output directory exists
    minimap2 -ax map-ont "$output_dir/${base_name}.fa" "$file" > "$output_dir/${base_name}.sam"
done


## moving some clusters into a phage subset folder - now done in step 8
#cd RNAquarium_outputs/virus_outputs/clusters_forminimap
#mkdir phage2
#mv *phage* phage2/
#mv *Phage* phage2/
#mv *Caudo* phage2/
#mv *Microvirid* phage2/
#mv *Carnation_lat* phage2/
#mv *Porcine_bastr* phage2/
#mv *Porcine_picobirn* phage2/
#mv *Viruses* phage2/
#mv *MAG_TPA* phage2/


## finally can move some interesting subsets into subset folder
## mkdir subsets
## mv *flu* subsets/
# scp -r eric.waltari@login01.czbiohub.org://hpc/projects/balla_group/sra_experiments/RNAquarium_60k_nonhostpipeline/RNAquarium_outputs/virus_outputs/clusters_forminimap clusters_2025


## rinse and repeat
#mv subsets/* .
#mv *Zebrafish* subsets/


## new tack - use sets from poster code
## then also separately make a SARS subfolder & also any Picorna subfolder

                                    
# mkdir clusters_forminimap_fishassoc
# cd clusters_forminimap

# mv *flue* ../clusters_forminimap_fishassoc/
# mv *Cyprinid* ../clusters_forminimap_fishassoc/
# mv *Cyvirus* ../clusters_forminimap_fishassoc/
# mv *Cyprivirus* ../clusters_forminimap_fishassoc/
# mv *Zebrafish* ../clusters_forminimap_fishassoc/
# mv *Coleura* ../clusters_forminimap_fishassoc/
# mv *Barramundi* ../clusters_forminimap_fishassoc/
# mv *salmon* ../clusters_forminimap_fishassoc/
# mv *Myotis* ../clusters_forminimap_fishassoc/
# mv *WFRC1* ../clusters_forminimap_fishassoc/
# mv *halfbeak* ../clusters_forminimap_fishassoc/
# mv *Astrovirus_3* ../clusters_forminimap_fishassoc/
# mv *snakehead* ../clusters_forminimap_fishassoc/
# mv *Snakehead* ../clusters_forminimap_fishassoc/
# mv *shark* ../clusters_forminimap_fishassoc/
# mv *Cutthroat* ../clusters_forminimap_fishassoc/
# mv *ctenopharyngodontis* ../clusters_forminimap_fishassoc/
# mv *catfish* ../clusters_forminimap_fishassoc/
# mv *Carp* ../clusters_forminimap_fishassoc/
# mv *Ranavirus* ../clusters_forminimap_fishassoc/
# mv *cichlid* ../clusters_forminimap_fishassoc/
# mv *Chuvivirus* ../clusters_forminimap_fishassoc/
# mv *Hardyhead* ../clusters_forminimap_fishassoc/
# mv *Salarius* ../clusters_forminimap_fishassoc/
# mv *Fish-associated* ../clusters_forminimap_fishassoc/
# mv *Luscinia* ../clusters_forminimap_fishassoc/
# mv *Phoenicopterus* ../clusters_forminimap_fishassoc/
# mv *Psittacidae* ../clusters_forminimap_fishassoc/
# mv *Psittaciform* ../clusters_forminimap_fishassoc/
# mv *Sprivivirus* ../clusters_forminimap_fishassoc/
# mv *Pimephales* ../clusters_forminimap_fishassoc/
# mv *cavefish* ../clusters_forminimap_fishassoc/
# mv *croaker* ../clusters_forminimap_fishassoc/
# mv *anguilla1* ../clusters_forminimap_fishassoc/
# mv *Catfish* ../clusters_forminimap_fishassoc/
# mv *angelfish* ../clusters_forminimap_fishassoc/
# mv *turtle* ../clusters_forminimap_fishassoc/
# mv *Symphysodon* ../clusters_forminimap_fishassoc/
# mv *Isavirus* ../clusters_forminimap_fishassoc/
# mv *Salmon* ../clusters_forminimap_fishassoc/
# mv *spikefish* ../clusters_forminimap_fishassoc/
# mv *rainbowfish* ../clusters_forminimap_fishassoc/


# cd ..
# mkdir clusters_forminimap_sarscov2
  
# cd clusters_forminimap

# mv *respiratory* ../clusters_forminimap_sarscov2/


# mkdir clusters_forminimap_picornas
# cd clusters_forminimap
# mv *picorna* ../clusters_forminimap_picornas/
# mv *Picorna* ../clusters_forminimap_picornas/
# cd ..

# cd ..
# mkdir clusters_forminimap_other


# mv *uncultured* ../clusters_forminimap_other/
# ## last!
# mv *Viruses* ../clusters_forminimap_other/


# mv *human* ../clusters_forminimap_other/
# mv *Human* ../clusters_forminimap_other/
# mv *dna* ../clusters_forminimap_other/
# mv *DNA* ../clusters_forminimap_other/


# mv *lenti* ../clusters_forminimap_other/
# mv *Lenti* ../clusters_forminimap_other/
# mv *Mouse* ../clusters_forminimap_other/
# mv *Monkey* ../clusters_forminimap_other/
# mv *monkey* ../clusters_forminimap_other/

# ribovir
# mv *Ribovir* ../clusters_forminimap_other/
# mv *ribovir* ../clusters_forminimap_other/
# mv *Crocidura* ../clusters_forminimap_other/
# mv *Porc* ../clusters_forminimap_other/
# mv *coron* ../clusters_forminimap_other/
# mv *Adinto* ../clusters_forminimap_other/
# mv *adinto* ../clusters_forminimap_other/
# mv *endog* ../clusters_forminimap_other/
# mv *Endog* ../clusters_forminimap_other/


# mkdir clusters_forminimap_more
# cd clusters_forminimap
# mv *Pinus* ../clusters_forminimap_more/
# mv *Cato* ../clusters_forminimap_more/
# mv *latens* ../clusters_forminimap_more/
# mv *closter* ../clusters_forminimap_more/
# mv *Latens* ../clusters_forminimap_more/
# mv *Closter* ../clusters_forminimap_more/
# mv *cato* ../clusters_forminimap_more/
# mv *arla* ../clusters_forminimap_more/

# mkdir clusters_forminimap_moren
# cd clusters_forminimap
# mv *n* ../clusters_forminimap_moren/


# scp -r eric.waltari@login01.czbiohub.org://hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/metatranscriptome_all/RNAquarium_outputs/virus_outputs/clusters_forminimap_fishassoc clusters_oct2025_fishassoc
# #downloaded
# #moved

# scp -r eric.waltari@login01.czbiohub.org://hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/metatranscriptome_all/RNAquarium_outputs/virus_outputs/clusters_forminimap_more clusters_oct2025_more
# #downloaded
# #moved

# scp -r eric.waltari@login01.czbiohub.org://hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/metatranscriptome_all/RNAquarium_outputs/virus_outputs/clusters_forminimap clusters_oct2025_many
# #downloaded
# #moved

# scp -r eric.waltari@login01.czbiohub.org://hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/metatranscriptome_all/RNAquarium_outputs/virus_outputs/clusters_forminimap_moren clusters_oct2025_many2
# #downloaded
# #moved

# clusters_forminimap_other
# clusters_forminimap_picornas
# clusters_forminimap_sarscov2

# scp -r eric.waltari@login01.czbiohub.org://hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/metatranscriptome_all/RNAquarium_outputs/virus_outputs/clusters_forminimap_other clusters_oct2025_other
# #downloaded
# #moved
# scp -r eric.waltari@login01.czbiohub.org://hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/metatranscriptome_all/RNAquarium_outputs/virus_outputs/clusters_forminimap_sarscov2 clusters_oct2025_sarscov2
# #downloaded
# #moved

# scp -r eric.waltari@login01.czbiohub.org://hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/metatranscriptome_all/RNAquarium_outputs/virus_outputs/clusters_forminimap_picornas clusters_oct2025_picornas
# #downloaded
# #moved

