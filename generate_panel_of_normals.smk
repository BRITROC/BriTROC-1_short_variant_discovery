#!/bin/env python

configfile: 'config.yaml'

import pandas

metadata = pandas.read_table("germline_metadata.tsv").set_index("fk_barcode", drop=False)
barcodes = metadata.index.unique().tolist()

rule all:
	input: 
		expand('normal_vcfs/SLX-15676_{barcode}.vcf', barcode=barcodes), 
		expand('normal_vcfs/{barcode}_pileups.table', barcode=barcodes), 
		'normals_for_pon_vcf.args',
		'pon.vcf.gz'

rule mutect2_normal_only:
	input:
		reference_genome=config['reference_genome'],
		bam='../SLX/SLX-15676/bam/{barcode}.bwamem.bam',
		germline_resource='../sample_swaps/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf.bgz'
	output: 'normal_vcfs/SLX-15676_{barcode}.vcf'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk Mutect2 \
			--reference {input.reference_genome} \
			--input {input.bam} \
			--germline-resource {input.germline_resource} \
			--output {output}'

rule create_list_file:
	input: normal_vcfs=expand('normal_vcfs/SLX-15676_{barcode}.vcf', barcode=barcodes)
	output: 'normals_for_pon_vcf.args'
	shell: 'ls normal_vcfs/*.vcf > {output}' 

rule create_panel_of_normals:
	input: rules.create_list_file.output
	output: 'pon.vcf.gz'
	shell: '/home/bioinformatics/software/gatk/gatk-4.0.11.0/gatk CreateSomaticPanelOfNormals --vcfs {input} -O {output}'

rule normal_sample_pileup_summary:
	input:
		bam='../SLX/SLX-15676/bam/{barcode}.bwamem.bam',
		biallelic_snps='somatic-b37_small_exac_common_3.fix_chr_names.vcf'
	output: 'normal_vcfs/{barcode}_pileups.table'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk GetPileupSummaries \
			-I {input.bam} \
			-V {input.biallelic_snps} \
			-L {input.biallelic_snps} \
			-O {output}'
