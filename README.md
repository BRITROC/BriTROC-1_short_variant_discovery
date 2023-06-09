# BriTROC-1 variant calling

The code repository for variant discovery used in the Smith & Bradley et al. (2023) publication.

# Assumptions

The code repository assumes that all of the TAm-Seq fastq data for this project is stored in a directory called fastq within subdirectories named after the name of the respective pooled library, such that it obeys the following pattern:

'{pooled_libeary_id}/{pooled_library_id}.{index_id}.{flowcell_id}.'{lane_id}.r_{read_number}.fq.gz'

# dependencies

Most dependencies are managed via snakemake, so the code must be executed using snakemake installed within a conda environment

# execution

The workflow can be executed from the main directory using the command 'snakemake -j 1 -n'
