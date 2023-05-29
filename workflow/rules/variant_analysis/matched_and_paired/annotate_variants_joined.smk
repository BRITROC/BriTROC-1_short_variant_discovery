def get_relevant_patient_list(wildcards):
	if wildcards.analysis_type=='panel_6_28':
		return(matched_somatic_patients)
	elif wildcards.analysis_type=='panel_28_only':
		return(matched_somatic_patients_panel_28_only)

rule get_somatic_variants_only_from_matched_and_paired_analysis:
	input: rules.concat_vcfs.output
	output: 'results/variant_analysis/matched/{analysis_type}/{patient_id}.filtered2.somatics_only.vcf'
	shell: 'grep -e "SOMATIC" -e "#" {input} > {output}'

rule ensure_tech_rep_genotypes_match_matched_analysis:
	input: 
		combined_vcfs=rules.get_somatic_variants_only_from_matched_and_paired_analysis.output,
		germline_metadata='config/matched_germline_metadata.tsv',
		tumour_metadata='config/matched_somatic_metadata.tsv'
	params:
		variant_quality_score_threshold=500,
		C_to_G_maf_threshold=0.23,
		C_to_G_maf_diff_threshold=0.30,
		includes_germline_sample_column=True,
		includes_tumour_type_analysis=True,
		includes_germline_variants=False
	output: 
		archival_samples='results/variant_analysis/matched/{analysis_type}/paired/{patient_id}.archival.filtered3.vcf',
		relapse_samples='results/variant_analysis/matched/{analysis_type}/paired/{patient_id}.relapse.filtered3.vcf',
	script: '../../../scripts/annotate_variants_joined/view_square_vcfs.R'

rule collate_and_filter_tumour_type_specific_vcf_files:
	input: lambda wildcards: expand('results/variant_analysis/matched/{analysis_type}/paired/{patient_id}.{tumour_type}.filtered3.vcf', analysis_type=wildcards.analysis_type, patient_id=get_relevant_patient_list(wildcards), tumour_type=wildcards.tumour_type)
	output: 'results/variant_analysis/matched/{analysis_type}/paired/collated/{tumour_type}_filtered3_joined.tsv'
	script: '../../../scripts/annotate_variants_joined/filtered4_files_joined.R'

rule MTBP_filter_curated_paired_results:
	input: rules.collate_and_filter_tumour_type_specific_vcf_files.output
	output: 'results/variant_analysis/matched/{analysis_type}/paired/collated/filtered3_{tumour_type}_joined_MTBP_filtered.tsv'
	script: '../../../scripts/annotate_variants_joined/apply_MTBP_filter.R'

rule collate_and_filter_octopus_vep_files_matched_analysis:
	input: 
		vep_files= lambda wildcards: expand('results/variant_analysis/matched/{analysis_type}/{patient_id}.filtered.vep.vcf',patient_id=matched_somatic_patients, analysis_type=wildcards.analysis_type),
		vcf_file=rules.MTBP_filter_curated_paired_results.output
	output: 
		vep_output='results/variant_analysis/matched/{analysis_type}/paired/collated/filtered_vep_calls_octopus_joined_{tumour_type}.tsv',
		vep_reduced='results/variant_analysis/matched/{analysis_type}/paired/collated/BriTROC-1_matched_and_paired_somatic_variants_{tumour_type}.tsv',
		vcf_output='results/variant_analysis/matched/{analysis_type}/paired/collated/filtered_calls_octopus_joined_{tumour_type}.vcf'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'
