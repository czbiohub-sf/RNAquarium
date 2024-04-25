TMP_OUTPUT := $(shell mktemp -d)

.PHONY: test
test:
	nextflow run \
		-profile test,slurm,singularity \
		-name test_run \
		--accession_list test/accession_subset.txt \
		--bioproj_map false \
		--publish_dir $(TMP_OUTPUT) \
		main.nf
	@nextflow clean test_run -f
