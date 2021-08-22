rule all:
	input: 
		'config/germline_metadata.tsv',
		'config/somatic_metadata',
		'config/paired_somatic_metadata.tsv',
		'config/paired_somatic_panel_28_metadata.tsv'

# generate germline and somatic metadata for the examination of all genes found in both panel 6 and panel 28
rule get_metadata:
	input:
	output: 
		germline_metadata='config/germline_metadata.tsv',
		somatic_metadata='config/somatic_metadata.tsv'
	script: '../scripts/generate_metadata/database_joins.R'

# generate germline and somatic metadata for the exmaination of all genes found in both panel 6 and panel 28 - ensuring there is at least one sample of each type
rule get_matched_metadata:
	input:
	output: 
		germline_metadata='config/matched_germline_metadata.tsv',
		somatic_metadata='config/matched_somatic_metadata.tsv'
	script: '../scripts/generate_metadata/matched_database_joins.R'

# generate germline somatic metadata for the examination of all genes found in panel 28 but not panel 6 - ensuring there is at least one sample of each type
rule get_matched_metadata_panel_28_only:
	input:
	output: 
		germline_metadata='config/germline_metadata_panel_28_only.tsv',
		somatic_metadata='config/somatic_metadata_panel_28_only.tsv'
	script: '../scripts/generate_metadata/matched_database_joins_panel_28_only.R'

#######################

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

