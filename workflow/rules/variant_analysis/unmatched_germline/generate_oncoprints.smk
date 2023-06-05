rule curate_germline_variants:
	output: 
		filtered_germline_variants='results/variant_analysis/germline/ampliconseq_pipeline/{analysis_type}/collated/filtered_germline_variants.tsv',
		HGVS_output_file='results/variant_analysis/germline/ampliconseq_pipeline/{analysis_type}/collated/filtered_germline_variants_HGVS.txt', # for MTBP
		filtered_germline_variants_vcf='results/variant_analysis/germline/ampliconseq_pipeline/{analysis_type}/collated/filtered_germline_variants.vcf' # for MTBP
	params: gene_set_analysed=get_gene_set_analysed
	script: '../../../scripts/curate_germline_variants.R'

rule MTBP_filter_germline_variants:
	input: rules.curate_germline_variants.output.filtered_germline_variants
	params: gene_set_analysed=get_gene_set_analysed
	output: 'results/variant_analysis/germline/ampliconseq_pipeline/{analysis_type}/collated/final_germline_variants.tsv'
	script: '../../../scripts/MTBP_curate_germline_variants.R'

#rule generate_germline_oncoprints:
#	input: filtered_germline_variants=rules.curate_germline_variants.output.filtered_germline_variants
#	output: germline_oncoprint='plots/britroc_germline_oncoprint.pdf' 
#	script: '../scripts/generate_oncoprints/generate_germline_oncoprint.R'

rule generate_germline_oncoprints_octopus:
	input: 
		filtered_germline_variants=rules.collate_and_filter_octopus_vep_files_germline.output,
		DNA_sample_file_path='britroc1_db/database_text_file_output/sample.tsv',
		patient_table_file_path='britroc1_db/database_text_file_output/patient.tsv'
	output: germline_oncoprint='plots/{analysis_type}_britroc_germline_oncoprint_octopus.png' 
	script: '../../../scripts/generate_oncoprints/generate_germline_oncoprint_octopus.R'

rule generate_germline_oncoprints:
	input: 
		filtered_germline_variants='short_variants_SI.tsv',
		DNA_sample_file_path='britroc1_db/database_text_file_output/sample.tsv',
		patient_table_file_path='britroc1_db/database_text_file_output/patient.tsv'
	output: germline_oncoprint='plots/{analysis_type}_britroc_germline_oncoprint_paper_revision.png' 
	script: '../../../scripts/generate_oncoprints/generate_germline_oncoprint_octopus.R'
