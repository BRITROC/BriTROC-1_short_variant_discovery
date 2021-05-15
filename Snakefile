# Master snakefile for this snakemake repository

configfile: 'config/config.yaml'

import pandas
import os

germline_metadata = pandas.read_table("config/germline_metadata.tsv").set_index("fk_barcode", drop=False)
germline_barcodes = germline_metadata.index.unique().tolist()

germline_panel_28_metadata = germline_metadata[(germline_metadata.fk_amplicon_panel==28)]

paired_somatic_metadata = pandas.read_table("config/paired_somatic_metadata.tsv").set_index("fk_sample", drop=False)
paired_somatic_samples = paired_somatic_metadata.index.unique().tolist()

paired_somatic_panel_28_metadata = pandas.read_table("config/paired_somatic_panel_28_metadata.tsv").set_index("fk_sample", drop=False)
paired_somatic_panel_28_samples = paired_somatic_panel_28_metadata.index.unique().tolist()

paired_archival_samples = paired_somatic_metadata[(paired_somatic_metadata.type=='archival')]
paired_archival_samples = paired_archival_samples.index.unique().tolist()

paired_relapse_samples = paired_somatic_metadata[(paired_somatic_metadata.type=='relapse')]
paired_relapse_samples = paired_relapse_samples.index.unique().tolist()

all_somatic_metadata = pandas.read_table("config/somatic_metadata.tsv").set_index("fk_sample", drop=False)
all_somatic_samples = all_somatic_metadata.index.unique().tolist()

archival_metadata = all_somatic_metadata[(all_somatic_metadata.type=='archival')]
all_archival_samples = archival_metadata.index.unique().tolist()

relapse_metadata = all_somatic_metadata[(all_somatic_metadata.type=='relapse')]
all_relapse_samples = relapse_metadata.index.unique().tolist()

somatic_germline_28_metadata = pandas.read_table("config/somatic_samples_with_germline_28.tsv").set_index("fk_sample", drop=False)
somatic_germline_28_samples = somatic_germline_28_metadata.index.unique().tolist()

somatic_tp53_metadata = pandas.read_table("config/somatic_samples_tp53.tsv").set_index("fk_sample", drop=False)
somatic_tp53_samples = somatic_tp53_metadata.index.unique().tolist()

somatic_tp53_archival_metadata = somatic_tp53_metadata[(somatic_tp53_metadata.type=='archival')]  
somatic_tp53_archival_samples = somatic_tp53_archival_metadata.index.unique().tolist()

#include: 'workflow/rules/process_gnomad_data.smk'
include: 'workflow/rules/generate_panel_of_normals.smk'
#include: 'workflow/rules/rsync.smk'
include: 'workflow/rules/generate_intersected_amplicons.smk'
include: 'workflow/rules/clean_bams.smk'
#include: 'workflow/rules/mutect2.smk'
include: 'workflow/rules/vardict.smk'
include: 'workflow/rules/octopus_tp53.smk'
include: 'workflow/rules/filter_mutect2_calls.smk'
include: 'workflow/rules/annotate_variants.smk'

rule all:
	input:
		#expand('results/rsync/patient-{britroc_number}_{sample_type}_{sample_id}_{barcode}', zip, britroc_number=all_somatic_metadata['fk_britroc_number'],sample_type=all_somatic_metadata['type'], sample_id=all_somatic_metadata['fk_sample'], barcode=all_somatic_metadata['fk_barcode']),
		#expand('results/rsync/patient-{britroc_number}_{sample_type}_{sample_id}_{barcode}', zip, britroc_number=somatic_tp53_metadata['fk_britroc_number'],sample_type=somatic_tp53_metadata['type'], sample_id=somatic_tp53_metadata['fk_sample'], barcode=somatic_tp53_metadata['fk_barcode'])
		#expand('results/rsync/patient-{britroc_number}_{sample_type}_{sample_id}_{barcode}', zip, britroc_number=germline_metadata['fk_britroc_number'],sample_type=germline_metadata['type'], sample_id=germline_metadata['fk_sample'], barcode=germline_metadata['fk_barcode'])
		'results/filtered_archival_vep_calls_octopus.tsv',
		'results/filtered_relapse_vep_calls_octopus.tsv',
		'results/tp53_collated_MAFs.tsv'
