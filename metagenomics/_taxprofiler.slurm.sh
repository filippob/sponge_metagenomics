#!/bin/bash
#SBATCH --job-name=taxprofiler    # Job name
#SBATCH --ntasks=1                    # Run on a single task
#SBATCH --cpus-per-task=1             # Declare 1 CPUs per task
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --time=00:05:00               # Time limit hrs:min:sec
#SBATCH --output=taxprofiler_%j.log   # Standard output and error log

nextflow run nf-core/taxprofiler \
	-r 2.0.1 \
	--input samplesheet.csv \
	--databases databases.csv \
	--run_kraken2 \
	--run_bracken \
	--skip_preprocessing_qc \
	-profile singularity, test
	#--run_metaphlan \


