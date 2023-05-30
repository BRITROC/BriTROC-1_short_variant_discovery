# somatic variants now found at the following path: results/variant_analysis/unmatched/paired/submission_results/collated/

rule prepare_data_for_whole_oncoprint_generation_unmatched_and_paired:
	input:
		filtered_non_TP53_variants_archival='results/variant_analysis/unmatched/{analysis_type}/paired/collated/archival_filtered_vep_calls_octopus_joined.targeted.tsv',
		filtered_non_TP53_variants_relapse='results/variant_analysis/unmatched/{analysis_type}/paired/collated/relapse_filtered_vep_calls_octopus_joined.targeted.tsv',
		filtered_TP53_variants_with_MAFs='results/variant_analysis/TP53/collated/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv',
	output: data_for_somatic_oncoprint='results/data_for_whole_cohort_oncoprint_{analysis_type}_intercalated.targeted.tsv'
	script: '../scripts/generate_oncoprints/prepare_data_for_somatic_oncoprints_cohort_both_tumour_types_only.R'

rule generate_unmatched_and_paired_oncoprint:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_whole_oncoprint_generation_unmatched_and_paired.output.data_for_somatic_oncoprint
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: somatic_oncoprint='plots/whole_cohort_oncoprints_{analysis_type}_intercalated.targeted.png' 
	script: '../scripts/generate_oncoprints/generate_somatic_oncoprints.R'

rule generate_unmatched_and_paired_oncoprint_diff_labels:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_whole_oncoprint_generation_unmatched_and_paired.output.data_for_somatic_oncoprint
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: somatic_oncoprint='plots/whole_cohort_oncoprints_{analysis_type}_intercalated.targeted.diff_labels.png' 
	script: '../scripts/generate_oncoprints/generate_somatic_oncoprints_diff_labels.R'
