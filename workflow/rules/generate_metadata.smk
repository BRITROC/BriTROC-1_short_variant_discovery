rule all:
	input: 'germline_metadata.tsv', 'somatic_metadata', 'paired_somatic_metadata.tsv'

rule get_metadata:
	input:
	output: 
		germline_metadata='germline_metadata.tsv',
		somatic_metadata='somatic_metadata.tsv'
	script: '../scripts/database_joins.R'

rule get_paired_metadata: # for paired primary and relapse samples
	input:
	output: 
		somatic_metadata='paired_somatic_metadata.tsv'
	script: '../scripts/database_joins_paired_samples.R'
