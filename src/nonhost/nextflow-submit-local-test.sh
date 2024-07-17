#!/bin/bash

set -o verbose
PATH=$PATH:$PWD/bin nextflow run main.nf \
					-params-file params.local.yaml \
					-profile testing,mamba
