#!/bin/env python

configfile: 'config.yaml'

import pandas
import os

metadata = pandas.read_table("somatic_metadata.tsv").set_index("fk_sample", drop=False)
samples = metadata.index.unique().tolist()

germline_metadata = pandas.read_table("germline_metadata.tsv").set_index("fk_britroc_number", drop=False)

rule all:
	input: 
		expand('tumour_vcfs/{sample}.vcf', sample=samples),
		expand('tumour_vcfs/{SLX}_{barcode}_{flowcell}_{lane}_pileups.table', zip, SLX=metadata['fk_slx'], barcode=metadata['fk_barcode'], flowcell=metadata['flowcell'], lane=metadata['lane']),
		expand('tumour_vcfs/{SLX}_{barcode}_{flowcell}_{lane}_contamination.table', zip, SLX=metadata['fk_slx'], barcode=metadata['fk_barcode'], flowcell=metadata['flowcell'], lane=metadata['lane']),
		expand('tumour_vcfs/{SLX}_{barcode}_{flowcell}_{lane}_segmentation.table', zip, SLX=metadata['fk_slx'], barcode=metadata['fk_barcode'], flowcell=metadata['flowcell'], lane=metadata['lane']),
		expand('tumour_vcfs/{sample}.filtered.vcf', sample=samples),
		expand('tumour_vcfs/{sample}.filtered2.vcf', sample=samples),
		expand('tumour_vcfs/{sample}.filtered3.vcf', sample=samples),
		expand('tumour_vcfs/{sample}.filtered.vep.vcf', sample=samples)

def get_tumour_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	metadata_filt = metadata[(metadata.fk_sample == wildcards.sample)].head(n=2)

	SLX = metadata_filt.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	barcodes = metadata_filt.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()
	lane = metadata_filt.set_index('lane', drop=False)
	lane = lane.index.unique().tolist()
	flowcell = metadata_filt.set_index('fk_run', drop=False)
	flowcell = flowcell.index.unique().tolist()
	flowcell = flowcell[0].split('_')

	bam_file_1 = '../SLX/{}/bam/{}.{}.{}.s_{}.bam'.format(SLX[0],  ''.join(SLX), barcodes[0], flowcell[-1], lane[0])
	bam_file_2 = '../SLX/{}/bam/{}.{}.{}.s_{}.bam'.format(SLX[0],  ''.join(SLX), barcodes[1], flowcell[-1], lane[0])

	return([bam_file_1, bam_file_2])

def get_normal_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	metadata_filt = metadata[(metadata.fk_sample==wildcards.sample)].head(n=2)
	britroc_number = metadata_filt.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	germline_metadata_filt = germline_metadata[germline_metadata['fk_britroc_number'] == britroc_number[0]]

	SLX = germline_metadata_filt.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()

	barcodes = germline_metadata_filt.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()

	lane = germline_metadata_filt.set_index('lane', drop=False)
	lane = lane.index.unique().tolist()
	
	flowcell = germline_metadata_filt.set_index('fk_run', drop=False)
	flowcell = flowcell.index.unique().tolist()
	flowcell = flowcell[0].split('_')

	if SLX[0] == 'SLX-9629':

		bam_file_1 = '../SLX/{}/samplebam/{}.{}.bam'.format(SLX[0], SLX[0], barcodes[0])
		return([bam_file_1])

	else:
		bam_file_1 = '../SLX/{}/bam/{}.{}.{}.s_{}.bam'.format(SLX[0],  ''.join(SLX), barcodes[0], flowcell[-1], lane[0])
		bam_file_2 = '../SLX/{}/bam/{}.{}.{}.s_{}.bam'.format(SLX[0],  ''.join(SLX), barcodes[1], flowcell[-1], lane[0])

		return([bam_file_1, bam_file_2])

def get_normal_samples(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	metadata_filt = metadata[(metadata.fk_sample==wildcards.sample)].head(n=2)
	britroc_number = metadata_filt.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	germline_metadata_filt = germline_metadata[germline_metadata['fk_britroc_number'] == britroc_number[0]]

	samples = germline_metadata_filt.set_index('fk_sample', drop=False)
	samples = samples.index.unique().tolist()	

	SLX = germline_metadata_filt.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()

	if SLX[0] == 'SLX-9629':
		pass
		#samples[0] = samples[0].replace('-','')
	else:
		samples.append(samples[0] + '_d')
		# this is accounting for a bug in the formatting of BAM headers
		#samples[0] = samples[0].replace('-','')
		#samples[1] = samples[1].replace('-','')

	return(samples)

rule get_tumour_pileup_summary:
	input:
		bam='../SLX/{slx}/bam/{slx}.{barcode}.{flowcell}.s_{lane}.bam',
		biallelic_snps='somatic-b37_small_exac_common_3.fix_chr_names.vcf'
	output: 'tumour_vcfs/{slx}_{barcode}_{flowcell}_{lane}_pileups.table'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk GetPileupSummaries \
			-I {input.bam} \
			-V {input.biallelic_snps} \
			-L {input.biallelic_snps} \
			-O {output}'

rule get_contamination_table:
	input: rules.get_tumour_pileup_summary.output,
	output: 
		contamination_table='tumour_vcfs/{slx}_{barcode}_{flowcell}_{lane}_contamination.table',
		tumour_segmentation='tumour_vcfs/{slx}_{barcode}_{flowcell}_{lane}_segmentation.table'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk CalculateContamination \
			-I {input} \
			--tumor-segmentation {output.tumour_segmentation} \
			-O {output.contamination_table}'
rule mutect2:
	input:
		reference_genome=config['reference_genome'],
		panel_of_normals='pon.vcf.gz',
		germline_resource='../sample_swaps/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf.bgz',
		interval_file='intersected_panel_6_28_amplicons.interval_list',
		tumour_bams=get_tumour_bam_files,
		normal_bams=get_normal_bam_files
	output: 
		tumour_vcf='tumour_vcfs/{sample}.vcf',
		f1r2='tumour_vcfs/{sample}_f1r2.tar.gz',
		bam_output='tumour_vcfs/{sample}.bam'
	threads: 4
	params:
		normal_sample_identifier=get_normal_samples
	run:
		print(input.normal_bams)
		if len(params.normal_sample_identifier) == 2:
			print('2 normal sample bam files')
			shell ('/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk Mutect2 \
				--input {input.tumour_bams[0]} \
				--input {input.tumour_bams[1]} \
				--input {input.normal_bams[0]} \
				--input {input.normal_bams[1]} \
				--normal-sample {params.normal_sample_identifier[0]} \
				--normal-sample {params.normal_sample_identifier[1]} \
				--intervals {input.interval_file} \
				--panel-of-normals {input.panel_of_normals} \
				--germline-resource {input.germline_resource} \
				--output {output.tumour_vcf} \
				--f1r2-tar-gz {output.f1r2} \
				--minimum-allele-fraction 0.00 \
				--bam-output {output.bam_output} \
				--reference {input.reference_genome}')
		else:
			print('1 normal sample bam file')
			shell ('/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk Mutect2 \
				--input {input.tumour_bams[0]} \
				--input {input.tumour_bams[1]} \
				--input {input.normal_bams[0]} \
				--normal-sample {params.normal_sample_identifier[0]} \
				--intervals {input.interval_file} \
				--panel-of-normals {input.panel_of_normals} \
				--germline-resource {input.germline_resource} \
				--output {output.tumour_vcf} \
				--f1r2-tar-gz {output.f1r2} \
				--minimum-allele-fraction 0.00 \
				--bam-output {output.bam_output} \
				--reference {input.reference_genome}')

def get_contamination_tables(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	metadata_filt = metadata[(metadata.fk_sample==wildcards.sample)].head(n=2)
	SLX = metadata_filt.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	barcodes = metadata_filt.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()
	lane = metadata_filt.set_index('lane', drop=False)
	lane = lane.index.unique().tolist()
	flowcell = metadata_filt.set_index('fk_run', drop=False)
	flowcell = flowcell.index.unique().tolist()
	flowcell = flowcell[0].split('_')

	contamination_table_1 = 'tumour_vcfs/{}_{}_{}_{}_contamination.table'.format(SLX[0], barcodes[0], flowcell[-1], lane[0])
	contamination_table_2 = 'tumour_vcfs/{}_{}_{}_{}_contamination.table'.format(SLX[0], barcodes[1], flowcell[-1], lane[0])

	return([contamination_table_1, contamination_table_2])

def get_segmentation_tables(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	metadata_filt = metadata[(metadata.fk_sample==wildcards.sample)].head(n=2)
	SLX = metadata_filt.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	barcodes = metadata_filt.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()
	lane = metadata_filt.set_index('lane', drop=False)
	lane = lane.index.unique().tolist()
	flowcell = metadata_filt.set_index('fk_run', drop=False)
	flowcell = flowcell.index.unique().tolist()
	flowcell = flowcell[0].split('_')

	segmentation_table_1 = 'tumour_vcfs/{}_{}_{}_{}_segmentation.table'.format(SLX[0], barcodes[0], flowcell[-1], lane[0])
	segmentation_table_2 = 'tumour_vcfs/{}_{}_{}_{}_segmentation.table'.format(SLX[0], barcodes[1], flowcell[-1], lane[0])

	return([segmentation_table_1, segmentation_table_2])

rule LearnReadOrientationModel:
	input: 'tumour_vcfs/{sample}_f1r2.tar.gz'
	output: 'tumour_vcfs/{sample}_artifiact_priors.tar.gz'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk LearnReadOrientationModel \
			--input {input} \
			--output {output}'

rule mark_mutect_calls_for_filtering:
	input:
		reference=config['reference_genome'],
		mutect2_vcf=rules.mutect2.output.tumour_vcf,
		f1r2=rules.LearnReadOrientationModel.output,
		contamination_tables=get_contamination_tables,
		segmentation_tables=get_segmentation_tables
	output: 'tumour_vcfs/{sample}.filtered.vcf'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk FilterMutectCalls \
			-R {input.reference} \
			-V {input.mutect2_vcf} \
			--orientation-bias-artifact-priors {input.f1r2} \
			--contamination-table {input.contamination_tables[0]} \
			--contamination-table {input.contamination_tables[1]} \
			--tumor-segmentation {input.segmentation_tables[0]} \
			--tumor-segmentation {input.segmentation_tables[1]} \
			--min-allele-fraction 0.00 \
			-O {output}'

rule filter_variants:
	input: 
		filtered_vcf=rules.mark_mutect_calls_for_filtering.output,
		reference=config['reference_genome']
	output: 'tumour_vcfs/{sample}.filtered2.vcf'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk SelectVariants \
			-R {input.reference} \
			-V {input.filtered_vcf} \
			--exclude-filtered \
			-O {output}'

rule filter_variants2:
	input: 
		filtered_vcf=rules.filter_variants.output
	output: 'tumour_vcfs/{sample}.filtered3.vcf'
	shell: 'sed s/^chr//g {input.filtered_vcf} > {output}'

rule vep:
	input:
		filtered_vcf=rules.filter_variants2.output
	output: 'tumour_vcfs/{sample}.filtered.vep.vcf'
	shell: 'ensembl-vep/vep \
			-i {input.filtered_vcf} \
			-o {output} \
			--database \
			--force_overwrite \
			--check_existing \
			--sift b \
			--polyphen b \
			--canonical \
			--protein \
			--uniprot \
			--numbers \
			--domains \
			--variant_class \
			--biotype \
			--ccds \
			--symbol \
			--clin_sig_allele 0 \
			-a GRCh37 \
			--port 3337'


