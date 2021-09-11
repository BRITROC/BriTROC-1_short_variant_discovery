rule curate_germline_variants:
	output: filtered_germline_variants='results/final_germline_variants.tsv'
	script: '../scripts/curate_germline_variants.R'

#rule generate_germline_and_somatic_oncoprints:
#	input:
#		filtered_non_TP53_variants='results/filtered_{tumour_sample_type}_vep_calls_octopus_joined.tsv',
#		filtered_TP53_variants_with_MAFs='results/final_tp53/filtered_TP53_variants_with_MAFs.tsv',
#		clonality_status_of_TP53_variants='results/final_tp53/TP53_variants_with_clonality_classifications.tsv',
#		filtered_germline_variants=rules.curate_germline_variants.output.filtered_germline_variants
#	output: germline_and_somatic_oncoprint='plots/{tumour_sample_type}_germline_and_somatic_britroc_oncoprint.pdf' 
#	script: '../scripts/generate_oncoprints/generate_germline_and_somatic_oncoprint.R'

rule generate_somatic_oncoprints:
	input:
		filtered_non_TP53_variants='results/panel_28_only/filtered_{tumour_sample_type}_vep_calls_octopus_joined.tsv',
		filtered_TP53_variants_with_MAFs='results/final_tp53/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/final_tp53/TP53_variants_with_clonality_classifications.tsv',
	output: somatic_oncoprint='plots/{tumour_sample_type}_somatic_britroc_oncoprint_panel_28_only.pdf' 
	script: '../scripts/generate_oncoprints/generate_somatic_oncoprints_panel_28_only.R'

rule generate_somatic_oncoprints_intercalated:
	input:
		filtered_non_TP53_variants_archival='results/panel_28_only/filtered_archival_vep_calls_octopus_joined.tsv',
		filtered_non_TP53_variants_relapse='results/panel_28_only/filtered_relapse_vep_calls_octopus_joined.tsv',
		filtered_TP53_variants_with_MAFs='results/final_tp53/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/final_tp53/TP53_variants_with_clonality_classifications.tsv',
	output: somatic_oncoprint='plots/somatic_britroc_oncoprint_panel_28_only_intercalated.pdf' 
	script: '../scripts/generate_oncoprints/generate_somatic_oncoprints_panel_28_only_intercalated.R'

rule generate_somatic_oncoprints_intercalated2:
	input:
		filtered_non_TP53_variants_archival='results/filtered_archival_vep_calls_octopus_joined.tsv',
		filtered_non_TP53_variants_relapse='results/filtered_relapse_vep_calls_octopus_joined.tsv',
		filtered_TP53_variants_with_MAFs='results/final_tp53/filtered_TP53_variants_with_MAFs.tsv',
		clonality_status_of_TP53_variants='results/final_tp53/TP53_variants_with_clonality_classifications.tsv',
	output: somatic_oncoprint='plots/somatic_britroc_oncoprint_intercalated.pdf' 
	script: '../scripts/generate_oncoprints/generate_somatic_oncoprints_intercalated.R'

#rule generate_germline_oncoprints:
#	input:
#		filtered_germline_variants=rules.curate_germline_variants.output.filtered_germline_variants
#	output: germline_oncoprint='plots/britroc_germline_oncoprint.pdf' 
#	script: '../scripts/generate_oncoprints/generate_germline_oncoprint.R'
