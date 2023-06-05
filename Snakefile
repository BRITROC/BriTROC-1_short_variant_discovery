# Master snakefile for this snakemake repository

configfile: 'config/config.yaml'

import pandas
import os

# TP53 sequencing metadata
TP53_tumour_sequencing_metadata = pandas.read_table("config/TP53_tumour_amplicon_sequencing_metadata.tsv").set_index('fk_sample', drop=False)
TP53_sequenced_DNA_samples = TP53_tumour_sequencing_metadata.index.unique().tolist()

# nontumour sequencing metadata
nontumour_sequencing_metadata = pandas.read_table("config/nontumour_amplicon_sequencing_metadata.tsv").set_index("fk_britroc_number", drop=False)
patients_with_nontumour_sample_sequencing = nontumour_sequencing_metadata.index.unique().tolist()

# panel 28 sequencing metadata
# unmatched and unpaired analysis metadata
# panel 28 = TP53 + other important HGSC related genes
panel_28_tumour_sequencing_metadata = pandas.read_table("config/panel_28_tumour_amplicon_sequencing_metadata.tsv").set_index("fk_britroc_number", drop=False)
patients_with_panel_28_tumour_sequencing = panel_28_tumour_sequencing_metadata.index.unique().tolist()

# matched and paired analysis metadata
matched_and_paired_sequencing_metadata = pandas.read_table("config/matched_and_paired_somatic_metadata.tsv").set_index("fk_britroc_number", drop=False)
patients_with_matched_and_paired_sequencing = matched_and_paired_sequencing_metadata.index.unique().tolist()

# matched and unpaired analysis metadata
matched_and_unpaired_sequencing_metadata = pandas.read_table('config/somatic_metadata_panel_matched_and_unpaired.tsv').set_index("fk_britroc_number", drop=False)
patients_with_matched_and_unpaired_sequencing = matched_and_unpaired_sequencing_metadata.index.unique().tolist()

# unmatched and paired analysis metadata
unmatched_and_paired_sequencing_metadata = pandas.read_table('config/tumour_metadata_with_one_of_both_types.tsv').set_index("fk_britroc_number", drop=False)
patients_with_unmatched_and_paired_sequencing = unmatched_and_paired_sequencing_metadata.index.unique().tolist()

wildcard_constraints:
	nonoverlapping_id='[1-9]',
	patient_id='[0-9]{1,3}',
	sample='(JBLAB-[0-9]+|IM_[0-9]+)'

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
		'plots/whole_cohort_oncoprints_panel_28_intercalated.targeted.png'	
