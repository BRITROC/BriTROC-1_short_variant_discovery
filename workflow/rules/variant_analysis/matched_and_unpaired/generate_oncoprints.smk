def get_gene_set_analysed(wildcards):
	if wildcards.analysis_type == 'panel_6_28':
		return(['TP53','BRCA1','BRCA2','FANCM','BARD1','RAD51B','RAD51C','RAD51D','BRIP1','PALB2'])
	elif wildcards.analysis_type == 'panel_28_only':
		return(['TP53','NRAS','PIK3CA','CTNNB1','EGFR','BRAF','PTEN','KRAS','RB1','CDK12','NF1'])

rule prepare_data_for_somatic_oncoprint_generation:
	input:
		filtered_non_TP53_variants='results/variant_analysis/matched/{analysis_type}/collated/filtered_vep_calls_octopus_joined.tsv',
		filtered_TP53_variants_with_MAFs='results/variant_analysis/TP53/collated/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv',
	output: data_for_somatic_oncoprint='results/data_for_somatic_oncoprint_{analysis_type}_matched_and_unpaired.tsv'
	script: '../../../scripts/generate_oncoprints/prepare_data_for_somatic_oncoprints_cohort.R'

rule generate_matched_and_unpaired_oncoprint_ggplot2:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_somatic_oncoprint_generation.output.data_for_somatic_oncoprint
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: somatic_oncoprint='plots/matched_and_unpaired_oncoprint_ggplot2_{analysis_type}.png' 
	script: '../../../scripts/generate_oncoprints/generate_somatic_oncoprint_matched_and_unpaired.R'
