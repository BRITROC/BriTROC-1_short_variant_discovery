#!/bin/env python

configfile: 'config.yaml'

import pandas

metadata = pandas.read_table("germline_metadata.tsv").set_index("fk_barcode", drop=False)
barcodes = metadata.index.unique().tolist()

rule all:
	input: 
		expand('normal_vcfs/{SLX}_{barcode}_{flowcell}_{lane}.vcf', zip, SLX=metadata['fk_slx'], barcode=metadata['fk_barcode'], flowcell=metadata['flowcell'], lane=metadata['lane']), 
		expand('normal_vcfs/{SLX}_{barcode}_{flowcell}_{lane}_pileups.table', zip, SLX=metadata['fk_slx'], barcode=metadata['fk_barcode'], flowcell=metadata['flowcell'], lane=metadata['lane']), 
		'pon.vcf.gz'

rule mutect2_normal_only:
	input:
		reference_genome=config['reference_genome'],
		bam='../SLX/{slx}/bam/{slx}.{barcode}.{flowcell}.s_{lane}.bam',
		germline_resource='../sample_swaps/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf.bgz',
		interval_file='intersected_panel_6_28_amplicons.interval_list'
	output: 'normal_vcfs/{slx}_{barcode}_{flowcell}_{lane}.vcf'
	threads: 4
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk Mutect2 \
			--reference {input.reference_genome} \
			--input {input.bam} \
			--intervals {input.interval_file} \
			--germline-resource {input.germline_resource} \
			--output {output}'

rule create_list_file:
	input: expand('normal_vcfs/{SLX}_{barcode}_{flowcell}_{lane}.vcf', zip, SLX=metadata['fk_slx'], barcode=metadata['fk_barcode'], flowcell=metadata['flowcell'], lane=metadata['lane'])
	output: 'normals_for_pon_vcf.args'
	shell: 'ls {input} > {output}' 

rule create_panel_of_normals:
	input: rules.create_list_file.output
	output: 'pon.vcf.gz'
	threads: 4
	shell: '/home/bioinformatics/software/gatk/gatk-4.0.11.0/gatk CreateSomaticPanelOfNormals --vcfs {input} -O {output}'

rule normal_sample_pileup_summary:
	input:
		bam='../SLX/{slx}/bam/{slx}.{barcode}.{flowcell}.s_{lane}.bam',
		biallelic_snps='somatic-b37_small_exac_common_3.fix_chr_names.vcf'
	output: 'normal_vcfs/{slx}_{barcode}_{flowcell}_{lane}_pileups.table'
	threads: 4
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk GetPileupSummaries \
			-I {input.bam} \
			-V {input.biallelic_snps} \
			-L {input.biallelic_snps} \
			-O {output}'
