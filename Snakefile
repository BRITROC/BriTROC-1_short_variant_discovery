# Master snakefile for this snakemake repository

configfile: 'config/config.yaml'

import pandas
import os

# read-in metadata sheets
matched_germline_metadata = pandas.read_table("config/matched_germline_metadata.tsv").set_index("fk_barcode", drop=False)
matched_somatic_metadata = pandas.read_table("config/matched_somatic_metadata.tsv").set_index("fk_britroc_number", drop=False)
matched_somatic_metadata_panel_28_only = pandas.read_table("config/somatic_metadata_panel_28_only.tsv").set_index("fk_britroc_number", drop=False)
somatic_tp53_metadata = pandas.read_table("config/somatic_samples_tp53.tsv").set_index('fk_sample', drop=False)
all_tumour_metadata = pandas.read_table("config/all_tumour_metadata.tsv").set_index("fk_britroc_number", drop=False)
tumour_metadata_patients_with_both_types = pandas.read_table('config/tumour_metadata_with_one_of_both_types.tsv').set_index("fk_britroc_number", drop=False)

matched_germline_barcodes = matched_germline_metadata.index.unique().tolist()
matched_somatic_patients = matched_somatic_metadata.index.unique().tolist()
matched_somatic_patients_panel_28_only = matched_somatic_metadata_panel_28_only.index.unique().tolist()
somatic_tp53_samples = somatic_tp53_metadata.index.unique().tolist()
all_tumour_sample_patients = all_tumour_metadata.index.unique().tolist()
all_patients_with_tumour_samples_of_both_types = tumour_metadata_patients_with_both_types.index.unique().tolist()

#wildcard_constraints:
#	nonoverlapping_id='[1-9]',
#	patient_id='[0-9]{1,3}$'

# preprocessing
include: 'workflow/rules/bam_preprocessing/generate_intersected_amplicons.smk'
include: 'workflow/rules/bam_preprocessing/clean_bams.smk'

# TP53 variant analysis
include: 'workflow/rules/variant_analysis/TP53/octopus_TP53.smk'
include: 'workflow/rules/variant_analysis/TP53/filter_mutect2_calls_TP53.smk'
include: 'workflow/rules/variant_analysis/TP53/annotate_variants_TP53.smk'

# non-TP53 variant analysis
#include: 'workflow/rules/variant_analysis/non_TP53/octopus_joint.smk'
#include: 'workflow/rules/variant_analysis/non_TP53/filter_mutect2_calls_joint.smk'
#include: 'workflow/rules/variant_analysis/non_TP53/annotate_variants_joined.smk'

# cohort level analysis snakemake
include: 'workflow/rules/variant_analysis/cohort/octopus_joint_cohort.smk'
include: 'workflow/rules/variant_analysis/cohort/filter_octopus_calls.smk'
include: 'workflow/rules/variant_analysis/cohort/annotate_variants_joined.smk'

# oncoprint generation
include: 'workflow/rules/generate_oncoprints_cohort.smk'

rule all:
	input:
		#'plots/panel_6_28/somatic_oncoprint.png',
		#'plots/panel_6_28/germline_and_somatic_oncoprint.png',
		#'plots/panel_28_only/somatic_oncoprint.png',
		#'plots/panel_28_only/germline_and_somatic_oncoprint.png'
		#expand('results/variant_analysis/cohort/{patient_id}.filtered3.vcf', patient_id=all_tumour_sample_patients)
		#'results/variant_analysis/cohort/collated/filtered3_joined.tsv',
		#'results/variant_analysis/cohort/collated/filtered_vep_calls_octopus_joined.tsv',
		#'results/data_for_somatic_oncoprint.tsv',
		#'plots/whole_cohort_somatic_oncoprint.png',
		'plots/whole_cohort_oncoprints_intercalated.png'
