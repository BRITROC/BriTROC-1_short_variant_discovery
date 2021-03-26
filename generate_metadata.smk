rule all:
	input: 'germline_metadata.tsv', 'somatic_metadata.tsv'

rule get_metadata:
	input:
	output: 
		germline_metadata='germline_metadata.tsv',
		somatic_metadata='somatic_metadata.tsv'
	script: 'database_joins.R'

rule get_paired_metadata: # for paired primary and relapse samples
	input:
	output: 
		somatic_metadata='paired_somatic_metadata.tsv'
	script: 'database_joins_paired_samples.R'

rule get_metadata_tp53:
	input:
	output:
		germline_metadata='germline_metadata_tp53.tsv',
		somatic_metadata='somatic_metadata_tp53.tsv'
	script: 'database_joins_tp53.R'
