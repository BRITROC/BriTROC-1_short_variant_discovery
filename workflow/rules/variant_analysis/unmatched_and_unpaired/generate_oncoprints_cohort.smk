# copyright Thomas Bradley 2023 ('thomas.bradley@cruk.cam.ac.uk')

rule prepare_data_for_whole_oncoprint_generation_not_intercalated:
	input:
		filtered_non_TP53_variants='results/variant_analysis/unmatched/{analysis_type}/collated/filtered_vep_calls_octopus_joined.tsv',
		filtered_TP53_variants_with_MAFs='results/variant_analysis/TP53/collated/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv',
		germline_variants='BriTROC-1_germline_variants.tsv',
		DNA_sample_file_path='britroc1_db/database_text_file_output/sample.tsv',
		slx_library_file_path='britroc1_db/database_text_file_output/slx_library.tsv',
		experiments_file_path='britroc1_db/database_text_file_output/experiment.tsv'
	conda: '../../../../config/postgres.yaml'
	output: data_for_somatic_oncoprint='results/data_for_somatic_oncoprint_{analysis_type}.tsv'
	script: '../../../scripts/generate_oncoprints/prepare_data_for_somatic_oncoprints_cohort.R'

rule generate_whole_cohort_oncoprint_not_intercalated_ggplot2:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_whole_oncoprint_generation_not_intercalated.output.data_for_somatic_oncoprint,
		patient_table_file_path='britroc1_db/database_text_file_output/patients.tsv'
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: somatic_oncoprint='plots/whole_cohort_oncoprints_{analysis_type}_not_intercalated_ggplot2.png' 
	script: '../../../scripts/generate_oncoprints/generate_somatic_oncoprints_cohort_ggplot2.R'
