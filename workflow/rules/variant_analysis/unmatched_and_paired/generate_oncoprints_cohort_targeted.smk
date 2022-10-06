def get_gene_set_analysed(wildcards):
	return(["BRCA1","BRCA2","RAD51C","RAD51D","RAD51B","BRIP1","FANCM","PALB2","BARD1","CDK12","EGFR","PTEN","TP53","KRAS","BRAF","PIK3CA","CTNNB1","NF1","RB1","NRAS"])

rule prepare_data_for_whole_oncoprint_generation_intercalated_targeted:
	input:
		filtered_non_TP53_variants_archival='results/variant_analysis/unmatched/collated/filtered_archival_vep_calls_octopus_joined.targeted.tsv',
		filtered_non_TP53_variants_relapse='results/variant_analysis/unmatched/collated/filtered_relapse_vep_calls_octopus_joined.targeted.tsv',
		filtered_TP53_variants_with_MAFs='results/variant_analysis/TP53/collated/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv',
	output: data_for_somatic_oncoprint='results/data_for_whole_cohort_oncoprint_intercalated.targeted.tsv'
	script: '../scripts/generate_oncoprints/prepare_data_for_somatic_oncoprints_cohort_both_tumour_types_only.R'

rule generate_whole_cohort_oncoprint_intercalated_targeted:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_whole_oncoprint_generation_intercalated_targeted.output.data_for_somatic_oncoprint
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: somatic_oncoprint='plots/whole_cohort_oncoprints_intercalated.targeted.png' 
	script: '../scripts/generate_oncoprints/generate_somatic_oncoprints.R'

rule generate_whole_cohort_oncoprint_intercalated_targeted_diff_labels:
	input:
		data_for_somatic_oncoprint=rules.prepare_data_for_whole_oncoprint_generation_intercalated_targeted.output.data_for_somatic_oncoprint
	params: 
		gene_set_analysed=get_gene_set_analysed
	output: somatic_oncoprint='plots/whole_cohort_oncoprints_intercalated.targeted.diff_labels.png' 
	script: '../scripts/generate_oncoprints/generate_somatic_oncoprints_diff_labels.R'
