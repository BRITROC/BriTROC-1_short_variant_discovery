# Master snakefile for this snakemake repository

configfile: 'config/config.yaml'

import pandas
import os

matched_germline_metadata = pandas.read_table("config/matched_germline_metadata.tsv").set_index("fk_barcode", drop=False)
matched_germline_barcodes = matched_germline_metadata.index.unique().tolist()

matched_somatic_metadata = pandas.read_table("config/matched_somatic_metadata.tsv").set_index("fk_sample", drop=False)
matched_somatic_samples = matched_somatic_metadata.index.unique().tolist()

matched_somatic_metadata = pandas.read_table("config/matched_somatic_metadata.tsv").set_index("fk_britroc_number", drop=False)
matched_somatic_patients = matched_somatic_metadata.index.unique().tolist()

matched_archival_metadata= matched_somatic_metadata[(matched_somatic_metadata.type=='archival')]
matched_archival_samples = matched_archival_metadata.index.unique().tolist()

include: 'workflow/rules/generate_intersected_amplicons.smk'
include: 'workflow/rules/clean_bams.smk'
#include: 'workflow/rules/octopus.smk'
include: 'workflow/rules/octopus_joint.smk'
include: 'workflow/rules/filter_mutect2_calls_joint.smk'
include: 'workflow/rules/annotate_variants_joined.smk'

rule all:
	input:
		#expand('results/tumour_sample_vcfs_octopus/{patient_id}.filtered3.vcf', patient_id=matched_somatic_patients),
		#expand('results/tumour_sample_vcfs_octopus/{patient_id}.library_MAFs.vcf', patient_id=matched_somatic_patients),
		#'results/archival_filtered3_joined.tsv',
		#'results/relapse_filtered3_joined.tsv',
		#'results/filtered3_joined.tsv',
		'results/filtered_archival_vep_calls_octopus_joined.tsv',
		'results/filtered_relapse_vep_calls_octopus_joined.tsv'
		#'results/archival_variants_with_two_reps.tsv'
		#expand('results/artifacts/{sample}.tsv', sample=matched_somatic_samples)
