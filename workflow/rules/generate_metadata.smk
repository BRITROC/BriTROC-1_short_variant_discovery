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
	script: '../scripts/database_joins.R'

rule get_paired_metadata: # for paired primary and relapse samples
	input:
	output: 
		somatic_metadata='config/paired_somatic_metadata.tsv'
	script: '../scripts/database_joins_paired_samples.R'

rule get_paired_panel_28metadata: # for paired primary and relapse samples
	input:
	output: 
		somatic_metadata='config/paired_somatic_panel_28_metadata.tsv'
	script: '../scripts/database_joins_paired_samples_panel_28.R'
