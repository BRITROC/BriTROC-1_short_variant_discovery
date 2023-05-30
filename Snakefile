# Master snakefile for this snakemake repository

configfile: 'config/config.yaml'

import pandas
import os

# read-in metadata sheets
matched_germline_metadata = pandas.read_table("config/matched_germline_metadata.tsv").set_index("fk_barcode", drop=False)
matched_somatic_metadata = pandas.read_table("config/matched_somatic_metadata.tsv").set_index("fk_britroc_number", drop=False)
matched_somatic_metadata_panel_28_only = pandas.read_table("config/somatic_metadata_panel_28_only.tsv").set_index("fk_britroc_number", drop=False)
germline_metadata = pandas.read_table("config/germline_metadata.tsv").set_index("fk_britroc_number", drop=False)

somatic_tp53_metadata = pandas.read_table("config/somatic_samples_tp53.tsv").set_index('fk_sample', drop=False)
#somatic_tp53_metadata_archival = somatic_tp53_metadata[(somatic_tp53_metadata.type == 'archival')]
#somatic_tp53_metadata_relapse = somatic_tp53_metadata[(somatic_tp53_metadata.type == 'relapse')]
all_tumour_metadata = pandas.read_table("config/all_tumour_metadata.tsv").set_index("fk_britroc_number", drop=False)
tumour_metadata_patients_with_both_types = pandas.read_table('config/tumour_metadata_with_one_of_both_types.tsv').set_index("fk_britroc_number", drop=False)

matched_and_unpaired_germline_metadata = pandas.read_table("config/germline_metadata_panel_matched_and_unpaired.tsv").set_index("fk_barcode", drop=False)
matched_and_unpaired_somatic_metadata = pandas.read_table('config/somatic_metadata_panel_matched_and_unpaired.tsv').set_index("fk_britroc_number", drop=False)

germline_patients = germline_metadata.index.unique().tolist()
matched_germline_barcodes = matched_germline_metadata.index.unique().tolist()
matched_somatic_patients = matched_somatic_metadata.index.unique().tolist()
matched_somatic_patients_panel_28_only = matched_somatic_metadata_panel_28_only.index.unique().tolist()
somatic_tp53_samples = somatic_tp53_metadata.index.unique().tolist()
all_tumour_sample_patients = all_tumour_metadata.index.unique().tolist()
all_patients_with_tumour_samples_of_both_types = tumour_metadata_patients_with_both_types.index.unique().tolist()
matched_and_unpaired_somatic_metadata_patients = matched_and_unpaired_somatic_metadata.index.unique().tolist()

print(matched_and_unpaired_somatic_metadata_patients)

#wildcard_constraints:
#	nonoverlapping_id='[1-9]',
#	patient_id='[0-9]{1,3}$'

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
#include: 'workflow/rules/variant_analysis/unmatched_germline/deep_variant.smk'
include: 'workflow/rules/variant_analysis/unmatched_germline/filter_octopus_calls.smk'
include: 'workflow/rules/variant_analysis/unmatched_germline/annotate_variants.smk'
include: 'workflow/rules/variant_analysis/unmatched_germline/generate_oncoprints.smk'

# matched and unpaired variant analysis
include: 'workflow/rules/variant_analysis/matched_and_unpaired/octopus_joint.smk'
include: 'workflow/rules/variant_analysis/matched_and_unpaired/filter_mutect2_calls_joint.smk'
#include: 'workflow/rules/variant_analysis/matched_and_unpaired/get_matched_germline_variants.smk'
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
		#'plots/panel_6_28/somatic_oncoprint.png',
		'results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv', # TP53 analysis output
		'results/variant_analysis/germline/ampliconseq_pipeline/panel_6_28/collated/final_germline_variants.tsv', # ampliconseq germline variants
		'plots/panel_6_28_britroc_germline_oncoprint_octopus.png', # octopus germline
		'plots/matched_and_unpaired_oncoprint_ggplot2_panel_6_28.png',
		'plots/panel_6_28/matched_and_paired_oncoprint_somatic_variants_only.png',
		'plots/panel_6_28/matched_and_paired_oncoprint_germline_and_somatic_variants.png',
		'plots/whole_cohort_oncoprints_panel_28_not_intercalated_ggplot2.png',
		'plots/whole_cohort_oncoprints_intercalated.targeted.png'
		#'plots/whole_cohort_oncoprints_not_intercalated_ggplot2.png'
		#'plots/matched_and_unpaired_oncoprint_ggplot2_panel_6_28.png'
		#'plots/panel_6_28/germline_and_somatic_oncoprint.png',
		#'plots/panel_28_only/somatic_oncoprint.png',
		#'plots/panel_28_only/germline_and_somatic_oncoprint.png'
		#'results/variant_analysis/germline/octopus_unmatched/panel_6_28/merged/collated/filtered_vep_calls_octopus_joined.tsv'
		#'results/variant_analysis/germline/octopus_unmatched/panel_6_28/merged/collated/filtered_vep_calls_octopus_joined.tsv'
		#expand('results/variant_analysis/germline/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.1.panel_6_28.passed.vcf', zip, SLX_ID=matched_germline_metadata['fk_slx'], barcodes=matched_germline_metadata['fk_barcode'], flowcell=matched_germline_metadata['flowcell'], lane=matched_germline_metadata['lane']),
		#expand('results/variant_analysis/germline/jointly_called.group_{group_number}.vcf', group_number=list(range(1,21,1))),
		#'results/variant_analysis/germline/jointly_called.vcf'
		#expand('results/variant_analysis/germline/octopus_unmatched/panel_6_28/merged/{patient_id}.vcf', patient_id=germline_patients)		
		#'results/variant_analysis/germline/octopus_unmatched/panel_6_28/merged/collated/filtered_vep_calls_octopus_joined.tsv',
		#'plots/panel_6_28_britroc_germline_oncoprint_octopus.png'
		#expand('results/variant_analysis/unmatched_germline/panel_6_28/deepvariant/{patient_id}.1.vcf.gz', patient_id=germline_patients),
		#expand('results/variant_analysis/unmatched_germline/panel_6_28/deepvariant/all_patients.{nonoverlapping_id}.g.vcf.gz', nonoverlapping_id=[1,2,3,4]),
		#expand('results/variant_analysis/unmatched_germline/panel_6_28/deepvariant/all_patients.g.vcf.gz'),
		#'results/variant_analysis/unmatched_germline/panel_6_28/deepvariant/all_patients.vep.vcf',
		#'results/variant_analysis/unmatched_germline/panel_6_28/deepvariant/all_patients.vep.filtered.vcf',
		#'results/variant_analysis/germline/ampliconseq_pipeline/panel_6_28/collated/final_germline_variants.tsv',
		#'results/variant_analysis/germline/octopus_unmatched/panel_6_28/merged/collated/filtered_vep_calls_octopus_joined.tsv'		
		#'results/variant_analysis/unmatched/collated/filtered3_joined_MTBP_filtered.tsv',
		#'results/variant_analysis/unmatched/collated/filtered_vep_calls_octopus_joined.tsv',
		#'results/data_for_somatic_oncoprint.tsv',
		#'plots/whole_cohort_oncoprints_not_intercalated_ggplot2.png',
		#'results/variant_analysis/matched/panel_6_28/paired/collated/filtered3_archival_joined_MTBP_filtered.tsv',
		#'results/variant_analysis/matched/panel_6_28/paired/collated/filtered3_relapse_joined_MTBP_filtered.tsv',
		#'results/variant_analysis/matched/panel_6_28/paired/collated/filtered_vep_calls_octopus_joined_archival.tsv',
		#'results/variant_analysis/matched/panel_6_28/paired/collated/filtered_vep_calls_octopus_joined_relapse.tsv',
		#'results/variant_analysis/unmatched/paired/collated/archival_filtered3_joined_MTBP_filtered.tsv',
		#'results/variant_analysis/unmatched/paired/collated/relapse_filtered3_joined_MTBP_filtered.tsv',
		#'results/variant_analysis/unmatched/paired/collated/archival_BriTROC-1_unmatched_and_paired_variants.tsv',
		#'results/variant_analysis/unmatched/paired/collated/relapse_BriTROC-1_unmatched_and_paired_variants.tsv',
		#expand('results/variant_analysis/unmatched/paired/{patient_id}.{nonoverlapping_id}.targeted.vcf', patient_id=all_patients_with_tumour_samples_of_both_types, nonoverlapping_id=[1,2,3,4]),
		#expand('results/variant_analysis/unmatched/paired/{patient_id}.filtered2.targeted.vcf', patient_id=all_patients_with_tumour_samples_of_both_types),
		#'results/variant_analysis/unmatched/paired/collated/archival_BriTROC-1_unmatched_and_paired_variants.targeted.tsv',
		#'results/variant_analysis/unmatched/paired/collated/relapse_BriTROC-1_unmatched_and_paired_variants.targeted.tsv'
		#expand('results/variant_analysis/unmatched/{patient_id}.filtered.vep.vcf',patient_id=all_tumour_sample_patients)
		#'results/variant_analysis/matched/panel_28_only/collated/filtered3_joined.tsv'

		#expand('results/variant_analysis/matched_and_unpaired/panel_6_28/{patient_id}.filtered2.vcf', patient_id=matched_and_unpaired_somatic_metadata_patients),
		#expand('results/variant_analysis/matched/panel_6_28/{patient_id}.interval_file.txt', patient_id=matched_somatic_patients),
#		expand('results/variant_analysis/matched/panel_6_28/paired/{patient_id}.filtered2.targeted.vcf', patient_id=matched_somatic_patients)
		#expand('results/variant_analysis/germline/octopus_matched_analysis/panel_6_28/{patient_id}.germline.vep_filtered.vcf', patient_id=matched_and_unpaired_somatic_metadata_patients)
		#'results/variant_analysis/germline/octopus_matched_analysis/panel_6_28/collated/filtered3_joined.tsv'
		#'results/variant_analysis/germline/octopus_matched_analysis/panel_6_28/collated/germline.vep_filtered.vcf',
		#'results/variant_analysis/matched/panel_6_28/collated/filtered3_joined.tsv',
		#'results/data_for_somatic_oncoprint_panel_6_28_matched_and_unpaired.tsv',
		#'plots/panel_6_28/matched_and_unpaired_germline_and_somatic_oncoprint_matched_octopus_germline.png',
		#'plots/panel_6_28/matched_and_unpaired_germline_and_somatic_oncoprint.png',
		#'config/germline_metadata.tsv'
		#'results/variant_analysis/matched_and_unpaired/panel_6_28/collated/filtered_vep_calls_octopus_joined.tsv',
		#'results/variant_analysis/matched_and_unpaired/panel_28_only/collated/filtered_vep_calls_octopus_joined.tsv'
		#expand('results/variant_analysis/unmatched/{patient_id}.{nono.targeted.vcf', patient_id=all_patients_with_tumour_samples_of_both_types, nonoverlapping_id=[1,2,3,4])
		#'results/variant_analysis/unmatched/collated/filtered3_joined.tsv',
		#'results/variant_analysis/unmatched/collated/filtered_vep_calls_octopus_joined.tsv',
		#'results/data_for_somatic_oncoprint.tsv',
		#'plots/whole_cohort_somatic_oncoprint.png',
		#'plots/whole_cohort_oncoprints_intercalated.png',
		#'plots/whole_cohort_oncoprints_intercalated.targeted.png'
		#'plots/matched_and_unpaired_oncoprint_ggplot2_panel_6_28.png'
