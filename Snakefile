#!/bin/env python

configfile: 'config/config.yaml'

import pandas
import os

germline_metadata = pandas.read_table("config/germline_metadata.tsv").set_index("fk_barcode", drop=False)
germline_barcodes = germline_metadata.index.unique().tolist()

germline_panel_28_metadata = germline_metadata[(germline_metadata.fk_amplicon_panel==28)]

paired_somatic_metadata = pandas.read_table("config/paired_somatic_metadata.tsv").set_index("fk_sample", drop=False)
paired_somatic_samples = paired_somatic_metadata.index.unique().tolist()

paired_archival_samples = paired_somatic_metadata[(paired_somatic_metadata.type=='archival')]
paired_archival_samples = paired_archival_samples.index.unique().tolist()

paired_relapse_samples = paired_somatic_metadata[(paired_somatic_metadata.type=='relapse')]
paired_relapse_samples = paired_relapse_samples.index.unique().tolist()

include: 'workflow/rules/process_gnomad_data.smk'
include: 'workflow/rules/generate_panel_of_normals.smk'
include: 'workflow/rules/mutect2.smk'
include: 'workflow/rules/filter_mutect2_calls.smk'
include: 'workflow/rules/annotate_variants.smk'

rule all:
	input: 
		expand('results/tumour_sample_vcfs/{sample}.filtered.vep.vcf', sample=paired_somatic_samples),
		'results/filtered_archival_vep_calls.tsv',
		#'results/filtered_relapse_vep_calls.tsv',
		#'results/shared_calls.tsv'
