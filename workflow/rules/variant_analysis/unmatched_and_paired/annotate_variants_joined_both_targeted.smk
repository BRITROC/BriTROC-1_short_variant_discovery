rule ensure_tech_rep_genotypes_match_unmatched_targeted:
	input: 
		combined_vcfs=rules.concat_vcfs_octopus_unmatched_targeted.output,
		tumour_metadata='config/tumour_metadata_with_one_of_both_types.tsv'
	params:
		variant_quality_score_threshold=500,
		C_to_G_maf_threshold=0.23,
		C_to_G_maf_diff_threshold=0.30,
		includes_germline_sample_column=False,
		includes_tumour_type_analysis=True,
		includes_germline_variants=True
	output: 
		archival_samples='results/variant_analysis/unmatched/{analysis_type}/paired/{patient_id}.archival.filtered3.targeted.vcf',
		relapse_samples='results/variant_analysis/unmatched/{analysis_type}/paired/{patient_id}.relapse.filtered3.targeted.vcf',
	script: '../../../scripts/annotate_variants_joined/view_square_vcfs.R'

rule collate_unmatched_tumour_type_specific_vcf_files:
	input: lambda wildcards: expand('results/variant_analysis/unmatched/{analysis_type}/paired/{patient_id}.{tumour_type}.filtered3.targeted.vcf', patient_id=patients_with_panel_28_tumour_sequencing, tumour_type=wildcards.tumour_type, analysis_type=wildcards.analysis_type)
	output: 'results/variant_analysis/unmatched/{analysis_type}/paired/collated/{tumour_type}_filtered3_joined.targeted.tsv'
	script: '../../../scripts/annotate_variants_joined/filtered4_files_joined.R'

rule MTBP_filter_unmatched_tumour_type_specific_variants:
	input: rules.collate_unmatched_tumour_type_specific_vcf_files.output
	output: 'results/variant_analysis/unmatched/{analysis_type}/paired/collated/{tumour_type}_filtered3_joined_MTBP_filtered.targeted.tsv'
	script: '../../../scripts/annotate_variants_joined/apply_MTBP_filter_unmatched_and_unpaired.R'

rule collate_and_filter_octopus_vep_files_with_tumour_type_targeted:
	input:
		vep_files=lambda wildcards: expand('results/variant_analysis/unmatched/{analysis_type}/{patient_id}.filtered.vep.vcf',patient_id=patients_with_panel_28_tumour_sequencing, analysis_type=wildcards.analysis_type),
		vcf_file=rules.MTBP_filter_unmatched_tumour_type_specific_variants.output
	output:
		vep_output='results/variant_analysis/unmatched/{analysis_type}/paired/collated/{tumour_type}_filtered_vep_calls_octopus_joined.targeted.tsv',
		vep_reduced='results/variant_analysis/unmatched/{analysis_type}/paired/collated/{tumour_type}_BriTROC-1_unmatched_and_paired_variants.targeted.tsv',
		vcf_output='results/variant_analysis/unmatched/{analysis_type}/paired/collated/{tumour_type}_filtered_calls_octopus_joined.targeted.vcf'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'
