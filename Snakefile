# Master snakefile for this snakemake repository

configfile: 'config/config.yaml'

import pandas
import os

# read-in metadata sheets

# TP53 metadata
somatic_tp53_metadata = pandas.read_table("config/somatic_samples_tp53.tsv").set_index('fk_sample', drop=False)
somatic_tp53_samples = somatic_tp53_metadata.index.unique().tolist()

# germline metadata
germline_metadata = pandas.read_table("config/germline_metadata.tsv").set_index("fk_britroc_number", drop=False)
germline_patients = germline_metadata.index.unique().tolist()

# matched and paired analysis metadata
matched_germline_metadata = pandas.read_table("config/matched_germline_metadata.tsv").set_index("fk_barcode", drop=False)
matched_somatic_metadata = pandas.read_table("config/matched_somatic_metadata.tsv").set_index("fk_britroc_number", drop=False)

matched_germline_barcodes = matched_germline_metadata.index.unique().tolist()
matched_somatic_patients = matched_somatic_metadata.index.unique().tolist()

# matched and unpaired analysis metadata
matched_and_unpaired_germline_metadata = pandas.read_table("config/germline_metadata_panel_matched_and_unpaired.tsv").set_index("fk_barcode", drop=False)
matched_and_unpaired_somatic_metadata = pandas.read_table('config/somatic_metadata_panel_matched_and_unpaired.tsv').set_index("fk_britroc_number", drop=False)

matched_and_unpaired_somatic_metadata_patients = matched_and_unpaired_somatic_metadata.index.unique().tolist()

# unmatched and unpaired analysis metadata
all_tumour_metadata = pandas.read_table("config/all_tumour_metadata.tsv").set_index("fk_britroc_number", drop=False)
all_tumour_sample_patients = all_tumour_metadata.index.unique().tolist()

# unmatched and paired analysis metadata
tumour_metadata_patients_with_both_types = pandas.read_table('config/tumour_metadata_with_one_of_both_types.tsv').set_index("fk_britroc_number", drop=False)
all_patients_with_tumour_samples_of_both_types = tumour_metadata_patients_with_both_types.index.unique().tolist()

# read in list of python functions created for this workflow
exec(open('snakemake_rule_functions.py').read())

# preprocessing
include: 'workflow/rules/bam_preprocessing/generate_intersected_amplicons.smk'
include: 'workflow/rules/bam_preprocessing/clean_bams.smk'

# TP53 variant analysis
include: 'workflow/rules/variant_analysis/TP53/octopus_TP53.smk'
include: 'workflow/rules/variant_analysis/TP53/annotate_variants_TP53.smk'

# germline workflow
include: 'workflow/rules/variant_analysis/unmatched_germline/octopus.smk'
include: 'workflow/rules/variant_analysis/unmatched_germline/filter_octopus_calls.smk'
include: 'workflow/rules/variant_analysis/unmatched_germline/annotate_variants.smk'
include: 'workflow/rules/variant_analysis/unmatched_germline/generate_oncoprints.smk'

# matched and unpaired variant analysis
include: 'workflow/rules/variant_analysis/matched_and_unpaired/octopus_joint.smk'
include: 'workflow/rules/variant_analysis/matched_and_unpaired/filter_mutect2_calls_joint.smk'
include: 'workflow/rules/variant_analysis/matched_and_unpaired/annotate_variants_joined.smk'
include: 'workflow/rules/variant_analysis/matched_and_unpaired/generate_oncoprints.smk'

# matched and paired variant analysis
include: 'workflow/rules/variant_analysis/matched_and_paired/octopus_join_targeted.smk'
include: 'workflow/rules/variant_analysis/matched_and_paired/annotate_variants_joined.smk'
include: 'workflow/rules/variant_analysis/matched_and_paired/filter_octopus_calls_targeted.smk'
include: 'workflow/rules/variant_analysis/matched_and_paired/generate_oncoprints.smk'

# unmatched and unpaired analyses
include: 'workflow/rules/variant_analysis/unmatched_and_unpaired/octopus_joint_cohort.smk'
include: 'workflow/rules/variant_analysis/unmatched_and_unpaired/filter_octopus_calls.smk'
include: 'workflow/rules/variant_analysis/unmatched_and_unpaired/annotate_variants_joined.smk'
include: 'workflow/rules/variant_analysis/unmatched_and_unpaired/generate_oncoprints_cohort.smk'

# unmatched and paired analysis
include: 'workflow/rules/variant_analysis/unmatched_and_paired/octopus_joined_targeted.smk'
include: 'workflow/rules/variant_analysis/unmatched_and_paired/filter_octopus_calls_targeted.smk'
include: 'workflow/rules/variant_analysis/unmatched_and_paired/annotate_variants_joined_both_targeted.smk'
include: 'workflow/rules/variant_analysis/unmatched_and_paired/generate_oncoprints_cohort_targeted.smk'

rule all:
	input:
		'results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv', # TP53 analysis output
		'results/variant_analysis/germline/ampliconseq_pipeline/panel_6_28/collated/final_germline_variants.tsv', # ampliconseq germline variants
		'plots/panel_6_28_britroc_germline_oncoprint_octopus.png', # octopus germline
		'plots/matched_and_unpaired_oncoprint_ggplot2_panel_6_28.png',
		'plots/panel_6_28/matched_and_paired_oncoprint_somatic_variants_only.png',
		'plots/panel_6_28/matched_and_paired_oncoprint_germline_and_somatic_variants.png',
		'plots/whole_cohort_oncoprints_panel_28_not_intercalated_ggplot2.png',
		'plots/whole_cohort_oncoprints_intercalated.targeted.png'	
