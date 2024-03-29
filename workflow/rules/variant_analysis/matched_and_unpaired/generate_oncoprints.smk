# copyright Thomas Bradley 2023 ('thomas.bradley@cruk.cam.ac.uk')

rule prepare_data_for_somatic_oncoprint_generation:
	input:
		filtered_non_TP53_variants='results/variant_analysis/matched/{analysis_type}/collated/filtered_vep_calls_octopus_joined.tsv',
		filtered_TP53_variants_with_MAFs='results/variant_analysis/TP53/collated/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv',
		DNA_sample_file_path='britroc1_db/database_text_file_output/sample.tsv',
		slx_library_file_path='britroc1_db/database_text_file_output/slx_library.tsv',
		experiments_file_path='britroc1_db/database_text_file_output/experiment.tsv'
	output: data_for_somatic_oncoprint='results/data_for_somatic_oncoprint_{analysis_type}_matched_and_unpaired.tsv'
	script: '../../../scripts/generate_oncoprints/prepare_data_for_somatic_oncoprints_cohort.R'

rule generate_matched_and_unpaired_oncoprint_ggplot2_somatic_and_germline:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_somatic_oncoprint_generation.output.data_for_somatic_oncoprint,
		germline_data=rules.curate_germline_variants.output.filtered_germline_variants,
		patient_table_file_path='britroc1_db/database_text_file_output/patients.tsv'
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: germline_and_somatic_oncoprint='plots/{analysis_type}/matched_and_unpaired_germline_and_somatic_oncoprint.png' 
	script: '../../../scripts/generate_oncoprints/generate_germline_and_somatic_oncoprint_ggplot2_2_unpaired_hc.R'

rule generate_matched_and_unpaired_oncoprint_ggplot2_somatic_and_germline_matched_octopus_germline:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_somatic_oncoprint_generation.output.data_for_somatic_oncoprint,
		germline_data=rules.collate_and_filter_octopus_vep_files_germline.output,
		patient_table_file_path='britroc1_db/database_text_file_output/patients.tsv'
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: germline_and_somatic_oncoprint='plots/{analysis_type}/matched_and_unpaired_germline_and_somatic_oncoprint_matched_octopus_germline.png' 
	script: '../../../scripts/generate_oncoprints/generate_germline_and_somatic_oncoprint_ggplot2_2_unpaired.R'

rule generate_matched_and_unpaired_oncoprint_ggplot2_somatic_only:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_somatic_oncoprint_generation.output.data_for_somatic_oncoprint,
		patient_table_file_path='britroc1_db/database_text_file_output/patients.tsv'
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: somatic_oncoprint='plots/matched_and_unpaired_oncoprint_ggplot2_{analysis_type}.png' 
	script: '../../../scripts/generate_oncoprints/generate_somatic_oncoprint_matched_and_unpaired.R'
