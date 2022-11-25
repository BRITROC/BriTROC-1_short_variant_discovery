def get_gene_set_analysed(wildcards):
	if wildcards.analysis_type == 'panel_6_28':
		return(['TP53','BRCA1','BRCA2','FANCM','BARD1','RAD51B','RAD51C','RAD51D','BRIP1','PALB2'])
	elif wildcards.analysis_type == 'panel_28_only':
		return(['TP53','NRAS','PIK3CA','CTNNB1','EGFR','BRAF','PTEN','KRAS','RB1','CDK12','NF1'])

rule prepare_data_for_somatic_oncoprint_generation:
	input:
		filtered_matched_and_paired_variants_archival='results/variant_analysis/matched_and_paired/{analysis_type}/collated/filtered_archival_vep_calls_octopus_joined.tsv',
		filtered_matched_and_paired_variants_relapse='results/variant_analysis/matched_and_paired/{analysis_type}/collated/filtered_relapse_vep_calls_octopus_joined.tsv',
		filtered_TP53_variants_with_MAFs='results/variant_analysis/TP53/collated/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv',
	output: data_for_somatic_oncoprint='results/data_for_somatic_oncoprint_{analysis_type}.tsv'
	script: '../scripts/generate_oncoprints/prepare_data_for_somatic_oncoprints.R'

rule generate_germline_and_somatic_oncoprints:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_somatic_oncoprint_generation.output.data_for_somatic_oncoprint,
		germline_data=rules.curate_germline_variants.output.filtered_germline_variants
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: germline_and_somatic_oncoprint='plots/{analysis_type}/germline_and_somatic_oncoprint.png' 
	script: '../scripts/generate_oncoprints/generate_germline_and_somatic_oncoprint.R'

rule generate_somatic_oncoprints_ggplot2:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_somatic_oncoprint_generation.output.data_for_somatic_oncoprint,
		germline_data=rules.curate_germline_variants.output.filtered_germline_variants
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: germline_and_somatic_oncoprint='plots/{analysis_type}/germline_and_somatic_oncoprint_ggplot2.png' 
	script: '../scripts/generate_oncoprints/generate_germline_and_somatic_oncoprint_ggplot2.R'

rule generate_somatic_oncoprints_ggplot2_2:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_somatic_oncoprint_generation.output.data_for_somatic_oncoprint,
		germline_data=rules.curate_germline_variants.output.filtered_germline_variants
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: germline_and_somatic_oncoprint='plots/{analysis_type}/germline_and_somatic_oncoprint_ggplot2_2.png' 
	script: '../scripts/generate_oncoprints/generate_germline_and_somatic_oncoprint_ggplot2_2.R'

rule generate_somatic_oncoprints_intercalated:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_somatic_oncoprint_generation.output.data_for_somatic_oncoprint
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: somatic_oncoprint='plots/{analysis_type}/somatic_oncoprint.png' 
	script: '../scripts/generate_oncoprints/generate_somatic_oncoprints.R'