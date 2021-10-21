rule curate_germline_variants:
	output: filtered_germline_variants='results/final_germline_variants.tsv'
	script: '../scripts/curate_germline_variants.R'

rule generate_somatic_oncoprints_intercalated:
	input:
		filtered_non_TP53_variants_archival='results/panel_28_only/filtered_archival_vep_calls_octopus_joined.tsv',
		filtered_non_TP53_variants_relapse='results/panel_28_only/filtered_relapse_vep_calls_octopus_joined.tsv',
		filtered_TP53_variants_with_MAFs='results/final_tp53/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/final_tp53/TP53_variants_with_clonality_classifications.tsv',
	output: somatic_oncoprint='plots/somatic_britroc_oncoprint_panel_28_only_intercalated.pdf' 
	script: '../scripts/generate_oncoprints/generate_somatic_oncoprints_panel_28_only_intercalated.R'
