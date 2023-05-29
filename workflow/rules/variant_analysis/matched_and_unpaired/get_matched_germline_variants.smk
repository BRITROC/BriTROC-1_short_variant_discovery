rule get_octopus_predicted_matched_germline_calls:
	input: rules.concat_vcfs.output
	output: 'results/variant_analysis/germline/octopus_matched_analysis/{analysis_type}/{patient_id}.germline.vcf'
	shell: 'grep -v SOMATIC {input} > {output}'

rule remove_germline_variants_not_passing_QC:
	input: rules.get_octopus_predicted_matched_germline_calls.output
	output: 'results/variant_analysis/germline/octopus_matched_analysis/{analysis_type}/{patient_id}.germline.filtered.vcf'
	script: '../../../scripts/filter_matched_analysis_germline_variants.R'

rule collate_and_filter_germline_vcf_files:
	input: lambda wildcards: expand('results/variant_analysis/germline/octopus_matched_analysis/{analysis_type}/{patient_id}.germline.filtered.vcf', patient_id=matched_and_unpaired_somatic_metadata_patients, analysis_type=wildcards.analysis_type)
	output: 'results/variant_analysis/germline/octopus_matched_analysis/{analysis_type}/collated/filtered3_joined.tsv'
	script: '../../../scripts/annotate_variants_joined/filtered4_files_joined_germline.R'

# contains really inefficient join operations - must fix at some point
rule collate_and_filter_octopus_vep_files_matched_germline:
	input: 
		vep_files= lambda wildcards: expand('results/variant_analysis/matched/{analysis_type}/{patient_id}.filtered.vep.vcf',patient_id=matched_and_unpaired_somatic_metadata_patients, analysis_type=wildcards.analysis_type),
		vcf_file=rules.collate_and_filter_germline_vcf_files.output
	output: 'results/variant_analysis/germline/octopus_matched_analysis/{analysis_type}/collated/germline.vep_filtered.vcf'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined_germline.R'
