.PHONY: test
test:
	nextflow run \
		-profile test,slurm,singularity \
		-w /hpc/scratch/group.swe/rnaquarium/work \
		--accession_list test/accession_subset.txt \
		--bioproj_map false main.nf
