TMP_OUTPUT := $(shell mktemp -d)

.PHONY: test
test:
	nextflow run main.nf \
		-profile test,slurm,singularity \
		-resume test_run \
		--accession_list test/accession_subset.txt \
		--bioproj_map false \
		--publish_dir $(TMP_OUTPUT)

.PHONY: small
small:
	nextflow run \
		-profile test,slurm,singularity \
		-resume \
		-ansi-log false \
		--publish_dir output/small_run \
		main.nf
