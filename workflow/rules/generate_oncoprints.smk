def get_gene_set_analysed(wildcards):
	if wildcards['HRD_or_nonHRD'] == 'HRD':
		return(['TP53','BRCA1','BRCA2','FANCM','BARD1','RAD51B','RAD51C','RAD51D','BRIP1','PALB2'])
	elif wildcards['HRD_or_nonHRD'] == 'nonHRD':
		return(['TP53','NRAS','PIK3CA','CTNNB1','EGFR','BRAF','PTEN','KRAS','RB1','CDK12','NF1'])

def get_relevant_nonTP53_variants_file_archival(wildcards):
	if wildcards['HRD_or_nonHRD'] == 'HRD':
		return('results/filtered_archival_vep_calls_octopus_joined.tsv')
	elif wildcards['HRD_or_nonHRD'] == 'nonHRD':
		return('results/panel_28_only/filtered_archival_vep_calls_octopus_joined.tsv')

def get_relevant_nonTP53_variants_file_relapse(wildcards):
	if wildcards['HRD_or_nonHRD'] == 'HRD':
		return('results/filtered_relapse_vep_calls_octopus_joined.tsv')
	elif wildcards['HRD_or_nonHRD'] == 'nonHRD':
		return('results/panel_28_only/filtered_relapse_vep_calls_octopus_joined.tsv')

rule curate_germline_variants:
	output: filtered_germline_variants='results/final_germline_variants.tsv'
	script: '../scripts/curate_germline_variants.R'

rule generate_germline_oncoprints:
	input: filtered_germline_variants=rules.curate_germline_variants.output.filtered_germline_variants
	output: germline_oncoprint='plots/britroc_germline_oncoprint.png' 
	script: '../scripts/generate_oncoprints/generate_germline_oncoprint.R'

rule prepare_data_for_somatic_oncoprint_generation:
	input:
		filtered_non_TP53_variants_archival=get_relevant_nonTP53_variants_file_archival,
		filtered_non_TP53_variants_relapse=get_relevant_nonTP53_variants_file_relapse,
		filtered_TP53_variants_with_MAFs='results/final_tp53/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/final_tp53/TP53_variants_with_clonality_classifications.tsv',
	output: data_for_somatic_oncoprint=temp('results/data_for_somatic_oncoprint_{HRD_or_nonHRD}.tsv')
	script: '../scripts/generate_oncoprints/prepare_data_for_somatic_oncoprints.R'

rule generate_germline_and_somatic_oncoprints:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_somatic_oncoprint_generation.output.data_for_somatic_oncoprint,
		germline_data='results/final_germline_variants.tsv'
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: germline_and_somatic_oncoprint='plots/{HRD_or_nonHRD}_germline_and_somatic_oncoprint.png' 
	script: '../scripts/generate_oncoprints/generate_germline_and_somatic_oncoprint.R'

rule generate_somatic_oncoprints_intercalated:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_somatic_oncoprint_generation.output.data_for_somatic_oncoprint
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: somatic_oncoprint='plots/{HRD_or_nonHRD}_somatic_oncoprint.png' 
	script: '../scripts/generate_oncoprints/generate_somatic_oncoprints.R'
