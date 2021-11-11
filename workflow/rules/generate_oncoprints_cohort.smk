def get_gene_set_analysed(wildcards):
	return(["BRCA1","BRCA2","RAD51C","RAD51D","RAD51B","BRIP1","FANCM","PALB2","BARD1","CDK12","EGFR","PTEN","TP53","KRAS","BRAF","PIK3CA","CTNNB1","NF1","RB1","NRAS"])

rule prepare_data_for_somatic_oncoprint_generation:
	input:
		filtered_non_TP53_variants='results/variant_analysis/cohort/collated/filtered_vep_calls_octopus_joined.tsv',
		filtered_TP53_variants_with_MAFs='results/variant_analysis/TP53/collated/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv',
	output: data_for_somatic_oncoprint='results/data_for_somatic_oncoprint.tsv'
	script: '../scripts/generate_oncoprints/prepare_data_for_somatic_oncoprints_cohort.R'

rule generate_whole_cohort_oncoprint_not_intercalated:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_somatic_oncoprint_generation.output.data_for_somatic_oncoprint
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: somatic_oncoprint='plots/whole_cohort_oncoprints_not_intercalated.png' 
	script: '../scripts/generate_oncoprints/generate_somatic_oncoprints_cohort.R'
