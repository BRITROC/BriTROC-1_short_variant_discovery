#!/bin/env python

#rule download_germline_resource:
#	input:
#	output: 'resources/gnomad.exomes.r2.1.1.sites.vcf.bgz'
#	shell: 'wget -O {output} https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz'

#rule decompress_gnomad_file:
#	input: rules.download_germline_resource.output
#	output: temp('resources/gnomad.exomes.r2.1.1.sites.vcf')
#	shell: 'gunzip -c {input} > {output}'

#rule add_chr_prefix:
#	input: rules.decompress_gnomad_file.output
#	output: temp('resources/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf')
#	shell: "sed -r 's/ID=([0-9XY]+)/ID=chr\\1/g' {input} | sed '/^#/!s/^/chr/1' > {output}" 

#rule compress_file:
#	input: rules.add_chr_prefix.output
#	output: protected('resources/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf.bgz')
#	shell: '/home/bioinformatics/software/htslib/htslib-1.6/bin/bgzip < {input} > {output}' 

#rule index_bgz_file:
#	input: rules.compress_file.output
#	output: protected('resources/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf.bgz.tbi')
#	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk IndexFeatureFile -I {input}'
