rule all:
	input: 
		'config/TP53_tumour_amplicon_sequencing_metadata.tsv',
		'config/panel_28_tumour_amplicon_sequencing_metadata.tsv',
		'config/nontumour_amplicon_sequencing_metadata.tsv',
		'config/tumour_metadata_with_one_of_both_types.tsv',
		'config/matched_and_paired_somatic_metadata.tsv',
		'config/somatic_metadata_panel_matched_and_unpaired.tsv'

configfile: 'config/config.yaml'

# generate metadata tables
rule collate_TP53_tumour_amplicon_sequencing_metadata:
	output: 'config/TP53_tumour_amplicon_sequencing_metadata.tsv'
	script: '../scripts/generate_metadata/database_joins_somatic_tp53.R'

rule collate_panel_28_tumour_amplicon_sequencing_metadata:
	input: rules.collate_TP53_tumour_amplicon_sequencing_metadata.output
	output: 'config/panel_28_tumour_amplicon_sequencing_metadata.tsv'
	script: '../scripts/generate_metadata/database_joins.R'

rule collate_nontumour_amplicon_sequencing_metadata:
	output: 'config/nontumour_amplicon_sequencing_metadata.tsv'
	script: '../scripts/generate_metadata/germline_metadata.R'

# generate patient lists

# generate metadata for all patients with at least one archival and relapse sample but not necessarily a germline sample
rule get_tumour_metadata_one_each_type:
	input: rules.collate_TP53_tumour_amplicon_sequencing_metadata.output
	output: 'config/tumour_metadata_with_one_of_both_types.tsv',
	script: '../scripts/generate_metadata/database_joins_one_tumour_sample_of_each_type.R'

# generate germline and somatic metadata for the exmaination of all genes found in both panel 6 and panel 28 - ensuring there is at least one sample of each type
rule get_matched_and_paired_metadata:
	output: 'config/matched_and_paired_somatic_metadata.tsv'
	script: '../scripts/generate_metadata/matched_database_joins.R'

rule get_matched_and_unpaired_only:
	output: 'config/somatic_metadata_panel_matched_and_unpaired.tsv'
	script: '../scripts/generate_metadata/matched_database_joins_matched_and_unpaired.R'
