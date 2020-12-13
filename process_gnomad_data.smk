#!/bin/env python

configfile: 'config.yaml'

import pandas

metadata = pandas.read_table("germline_metadata.tsv").set_index("fk_barcode", drop=False)
barcodes = metadata.index.unique().tolist()

rule all:
	input: '../sample_swaps/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf.bgz.tbi'

rule decompress_gnomad_file:
	input: '../sample_swaps/gnomad.exomes.r2.1.1.sites.vcf.bgz'
	output: temp('../sample_swaps/gnomad.exomes.r2.1.1.sites.vcf')
	shell: 'gunzip -c {input} > {output}'

rule add_chr_prefix:
	input: rules.decompress_gnomad_file.output
	output: temp('../sample_swaps/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf')
	shell: "sed -r 's/ID=([0-9XY]+)/ID=chr\\1/g' {input} | sed '/^#/!s/^/chr/1' > {output}" 

rule compress_file:
	input: rules.add_chr_prefix.output
	output: '../sample_swaps/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf.bgz'
	shell: '/home/bioinformatics/software/htslib/htslib-1.6/bin/bgzip < {input} > {output}' 

rule index_bgz_file:
	input: rules.compress_file.output
	output: '../sample_swaps/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf.bgz.tbi'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk IndexFeatureFile -I {input}'

rule add_chr_prefix_to_biallelic_SNPs:
	input: 'somatic-b37_small_exac_common_3.vcf'
	output: 'somatic-b37_small_exac_common_3.fix_chr_names.vcf'
	shell: "sed -r 's/ID=([0-9XY]+)/ID=chr\\1/g' {input} | sed '/^#/!s/^/chr/1' > {output}"

rule index_biallelic_snps_file:
	input: rules.add_chr_prefix_to_biallelic_SNPs.output
	output: 'somatic-b37_small_exac_common_3.fix_chr_names.vcf.idx'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk IndexFeatureFile -I {input}'

