# Master snakefile for this snakemake repository

configfile: 'config/config.yaml'

import pandas
import os

matched_germline_metadata = pandas.read_table("config/matched_germline_metadata.tsv").set_index("fk_barcode", drop=False)
matched_somatic_metadata = pandas.read_table("config/matched_somatic_metadata.tsv").set_index("fk_sample", drop=False)
matched_somatic_metadata = pandas.read_table("config/matched_somatic_metadata.tsv").set_index("fk_britroc_number", drop=False)
somatic_tp53_metadata = pandas.read_table("config/somatic_samples_tp53.tsv").set_index('fk_sample', drop=False)

matched_germline_barcodes = matched_germline_metadata.index.unique().tolist()
matched_somatic_samples = matched_somatic_metadata.index.unique().tolist()
matched_somatic_patients = matched_somatic_metadata.index.unique().tolist()
somatic_tp53_samples = somatic_tp53_metadata.index.unique().tolist()

include: 'workflow/rules/generate_intersected_amplicons.smk'
include: 'workflow/rules/clean_bams.smk'
include: 'workflow/rules/octopus_tp53.smk'
include: 'workflow/rules/filter_mutect2_calls.smk'
include: 'workflow/rules/annotate_variants_TP53.smk'

if config['panel_28_only'] == True:
	include: 'workflow/rules/joined/panel_28_only/octopus_joint_panel_28_only.smk'
	include: 'workflow/rules/joined/panel_28_only/annotate_variants_joined_panel_28_only.smk',
	include: 'workflow/rules/joined/panel_28_only/generate_oncoprints_panel_28_only.smk'

	rule all:
		input:
			'plots/somatic_britroc_oncoprint_panel_28_only_intercalated.pdf'

else if config['panel_28_only'] == False:
	include: 'workflow/rules/octopus_joint.smk'
	include: 'workflow/rules/filter_mutect2_calls_joint.smk'
	include: 'workflow/rules/generate_oncoprints.smk'

	rule all:
		input:
			'plots/HRD_somatic_oncoprint.pdf',
			'plots/HRD_germline_and_somatic_oncoprint.pdf'
