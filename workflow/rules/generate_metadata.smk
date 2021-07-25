rule all:
	input: 
		'config/germline_metadata.tsv',
		'config/somatic_metadata',
		'config/paired_somatic_metadata.tsv',
		'config/paired_somatic_panel_28_metadata.tsv'

rule get_metadata:
	input:
	output: 
		germline_metadata='config/germline_metadata.tsv',
		somatic_metadata='config/somatic_metadata.tsv'
	script: '../scripts/generate_metadata/database_joins.R'

rule get_matched_metadata:
	input:
	output: 
		germline_metadata='config/matched_germline_metadata.tsv',
		somatic_metadata='config/matched_somatic_metadata.tsv'
	script: '../scripts/generate_metadata/matched_database_joins.R'

rule get_paired_metadata: # for paired primary and relapse samples
	input:
	output: 
		somatic_metadata='config/paired_somatic_metadata.tsv'
	script: '../scripts/generate_metadata/database_joins_paired_samples.R'

rule get_paired_panel_28metadata: # for paired primary and relapse samples
	input:
	output: 
		somatic_metadata='config/paired_somatic_panel_28_metadata.tsv'
	script: '../scripts/generate_metadata/database_joins_paired_samples_panel_28.R'

rule get_somatic_samples_with_germline_28: # for paired primary and relapse samples
	input:
	output: 
		somatic_metadata='config/somatic_samples_with_germline_28.tsv'
	script: '../scripts/generate_metadata/database_joins_paired_somatic_samples_with_germline_28.R'

rule get_somatic_samples_tp53:
	input:
	output:
		somatic_metadata='config/somatic_samples_tp53.tsv'
	script: '../scripts/generate_metadata/database_joins_somatic_tp53.R'

