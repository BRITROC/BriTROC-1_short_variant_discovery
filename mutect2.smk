#!/bin/env python

configfile: 'config.yaml'

import pandas

metadata = pandas.read_table("somatic_metadata.tsv").set_index("fk_sample", drop=False)
samples = metadata.index.unique().tolist()

germline_metadata = pandas.read_table("germline_metadata.tsv").set_index("fk_britroc_number", drop=False)

rule all:
	input: 
		expand('tumour_vcfs/{sample}.vcf', sample=samples),
		metadata['table_file'],
		metadata['contamination_table'],
		metadata['segmentation_table'],
		expand('tumour_vcfs/{sample}.filtered.vcf', sample=samples),
		expand('tumour_vcfs/{sample}.filtered.annotated.maf', sample=samples)

def get_tumour_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	metadata_filt = metadata.filter(like=wildcards.sample, axis=0).head(n=2)
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
	metadata_filt = metadata.filter(like=wildcards.sample, axis=0).head(n=2)
	britroc_number = metadata_filt.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	germline_metadata_filt = germline_metadata.filter(like=str(britroc_number[0]), axis=0)

	SLX = germline_metadata_filt.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()

	barcodes = germline_metadata_filt.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()	

	bam_file_1 = '../SLX/{}/bam/{}.bwamem.bam'.format(SLX[0], barcodes[0])
	bam_file_2 = '../SLX/{}/bam/{}.bwamem.bam'.format(SLX[0], barcodes[1])

	return([bam_file_1, bam_file_2])

def get_normal_barcodes(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	metadata_filt = metadata.filter(like=wildcards.sample, axis=0).head(n=2)
	britroc_number = metadata_filt.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	germline_metadata_filt = germline_metadata.filter(like=str(britroc_number[0]), axis=0)

	barcodes = germline_metadata_filt.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()	

	return(barcodes)

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
		tumour_bams=get_tumour_bam_files,
		normal_bams=get_normal_bam_files
	output: 
		tumour_vcf='tumour_vcfs/{sample}.vcf',
		f1r2='tumour_vcfs/{sample}_f1r2.tar.gz'
	params:
		normal_sample_barcode=get_normal_barcodes
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk Mutect2 \
			--input {input.tumour_bams[0]} \
			--input {input.tumour_bams[1]} \
			--input {input.normal_bams[0]} \
			--input {input.normal_bams[1]} \
			--normal-sample {params.normal_sample_barcode[0]} \
			--normal-sample {params.normal_sample_barcode[1]} \
			--panel-of-normals {input.panel_of_normals} \
			--germline-resource {input.germline_resource} \
			--output {output.tumour_vcf} \
			--f1r2-tar-gz {output.f1r2} \
			--reference {input.reference_genome}'

def get_contamination_tables(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	metadata_filt = metadata.filter(like=wildcards.sample, axis=0).head(n=2)
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
	metadata_filt = metadata.filter(like=wildcards.sample, axis=0).head(n=2)
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

rule filter_mutect_calls:
	input:
		reference=config['reference_genome'],
		mutect2_vcf=rules.mutect2.output.tumour_vcf,
		f1r2=rules.mutect2.output.f1r2,
		contamination_tables=get_contamination_tables,
		segmentation_tables=get_segmentation_tables
	output: 'tumour_vcfs/{sample}.filtered.vcf'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk FilterMutectCalls \
			-R {input.reference} \
			-V {input.mutect2_vcf} \
			--contamination-table {input.contamination_tables[0]} \
			--contamination-table {input.contamination_tables[1]} \
			--tumor-segmentation {input.segmentation_tables[0]} \
			--tumor-segmentation {input.segmentation_tables[1]} \
			-O {output}'

rule funcotator:
	input: 
		reference=config['reference_genome'],
		filtered_vcf=rules.filter_mutect_calls.output
	output: 'tumour_vcfs/{sample}.filtered.annotated.maf'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk Funcotator \
			-R {input.reference} \
			-V {input.filtered_vcf} \
			--output-file-format MAF \
			--ref-version hg19 \
			--data-sources-path funcotator_dataSources.v1.6.20190124s \
			--remove-filtered-variants \
			-O {output}'

