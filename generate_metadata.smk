rule all:
	input: 'germline_metadata.tsv', 'somatic_metadata.tsv'

rule get_metadata:
	input:
	output: 
		germline_metadata='germline_metadata.tsv',
		somatic_metadata='somatic_metadata.tsv'
	script: 'database_joins.R'
