# somatic variants now found at the following path: results/variant_analysis/unmatched/paired/submission_results/collated/

# copyright Thomas Bradley 2023 ('thomas.bradley@cruk.cam.ac.uk')

# copyright Thomas Bradley 2023 ('thomas.bradley@cruk.cam.ac.uk')

rule prepare_data_for_whole_oncoprint_generation_unmatched_and_paired:
	input:
		filtered_non_TP53_variants_archival='results/variant_analysis/unmatched/{analysis_type}/paired/collated/archival_filtered_vep_calls_octopus_joined.targeted.tsv',
		filtered_non_TP53_variants_relapse='results/variant_analysis/unmatched/{analysis_type}/paired/collated/relapse_filtered_vep_calls_octopus_joined.targeted.tsv',
		filtered_TP53_variants_with_MAFs='results/variant_analysis/TP53/collated/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv',
		DNA_sample_file_path='britroc1_db/database_text_file_output/sample.tsv',
		slx_library_file_path='britroc1_db/database_text_file_output/slx_library.tsv',
		experiments_file_path='britroc1_db/database_text_file_output/experiment.tsv'
	output: data_for_somatic_oncoprint='results/data_for_whole_cohort_oncoprint_{analysis_type}_intercalated.targeted.tsv'
	script: '../scripts/generate_oncoprints/prepare_data_for_somatic_oncoprints_cohort_both_tumour_types_only.R'

rule generate_unmatched_and_paired_oncoprint:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_whole_oncoprint_generation_unmatched_and_paired.output.data_for_somatic_oncoprint,
		patient_table_file_path='britroc1_db/database_text_file_output/patients.tsv',
		chemotherapy_lines_table_file_path='britroc1_db/database_text_file_output/chemotherapy_lines.tsv'
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: somatic_oncoprint='plots/unmatched_and_paired_oncoprint_{analysis_type}.png' 
	script: '../scripts/generate_oncoprints/generate_unmatched_and_unpaired_oncoprint.R'
