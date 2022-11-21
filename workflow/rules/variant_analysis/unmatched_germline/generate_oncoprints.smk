def get_gene_set_analysed(wildcards):
	if wildcards.analysis_type == 'panel_6_28':
		return(['TP53','BRCA1','BRCA2','FANCM','BARD1','RAD51B','RAD51C','RAD51D','BRIP1','PALB2'])
	elif wildcards.analysis_type == 'panel_28_only':
		return(['TP53','NRAS','PIK3CA','CTNNB1','EGFR','BRAF','PTEN','KRAS','RB1','CDK12','NF1'])

rule curate_germline_variants:
	output: filtered_germline_variants='results/variant_analysis/germline/{analysis_type}/collated/final_germline_variants.tsv'
	params: gene_set_analysed=get_gene_set_analysed
	script: '../../../scripts/curate_germline_variants.R'

rule generate_germline_oncoprints:
	input: filtered_germline_variants=rules.curate_germline_variants.output.filtered_germline_variants
	output: germline_oncoprint='plots/britroc_germline_oncoprint.pdf' 
	script: '../scripts/generate_oncoprints/generate_germline_oncoprint.R'
