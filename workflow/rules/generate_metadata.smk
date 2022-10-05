rule all:
	input: 
		'config/matched_germline_metadata.tsv',
		'config/matched_somatic_metadata.tsv',
		'config/somatic_metadata_panel_28_only.tsv',
		'config/somatic_samples_tp53.tsv',
		'config/all_tumour_metadata.tsv'

configfile: 'config/config.yaml'
rule get_somatic_samples_tp53:
	output:
		somatic_metadata='config/somatic_samples_tp53.tsv'
	script: '../scripts/generate_metadata/database_joins_somatic_tp53.R'

# generate metadata for all tumour samples
rule get_all_tumour_metadata:
	input: somatic_tp53_metadata=rules.get_somatic_samples_tp53.output
	output: somatic_metadata='config/all_tumour_metadata.tsv'
	script: '../scripts/generate_metadata/database_joins.R'

# generate metadata for all patients with at least one archival and relapse sample but not necessarily a germline sample
rule get_tumour_metadata_one_each_type:
	input: somatic_tp53_metadata=rules.get_somatic_samples_tp53.output
	output: somatic_metadata='config/tumour_metadata_with_one_of_both_types.tsv',
	script: '../scripts/generate_metadata/database_joins_one_tumour_sample_of_each_type.R'

# generate germline and somatic metadata for the exmaination of all genes found in both panel 6 and panel 28 - ensuring there is at least one sample of each type
rule get_matched_metadata:
	output: 
		germline_metadata='config/matched_germline_metadata.tsv',
		somatic_metadata='config/matched_somatic_metadata.tsv'
	script: '../scripts/generate_metadata/matched_database_joins.R'

# generate germline somatic metadata for the examination of all genes found in panel 28 but not panel 6 - ensuring there is at least one sample of each type
rule get_matched_metadata_panel_28_only:
	output: 
		germline_metadata='config/germline_metadata_panel_28_only.tsv',
		somatic_metadata='config/somatic_metadata_panel_28_only.tsv'
	script: '../scripts/generate_metadata/matched_database_joins_panel_28_only.R'

rule get_matched_and_unpaired_only:
	output: 
		germline_metadata='config/germline_metadata_panel_matched_and_unpaired.tsv',
		somatic_metadata='config/somatic_metadata_panel_matched_and_unpaired.tsv'
	script: '../scripts/generate_metadata/matched_database_joins_matched_and_unpaired.R'
