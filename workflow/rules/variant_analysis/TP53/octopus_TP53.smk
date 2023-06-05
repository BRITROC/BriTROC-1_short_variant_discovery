rule convert_bed6_to_oct_format_tp53:
	input:   'resources/union_of_tp53_amplicons.targets.{nonoverlapping_id}.bed'                                     
	output: 'resources/union_of_tp53.{nonoverlapping_id}.targets.oct'
	script: '../scripts/octopus_formatting/convert_bed6_to_octopus.R'

rule octopus_tp53:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.convert_bed6_to_oct_format_tp53.output,
		tumour_bams=lambda wildcards: get_bam_files(wildcards, 'bam', 'tumour_TP53'),
		tumour_bam_indexes = lambda wildcards: get_bam_files(wildcards, 'bai', 'tumour_TP53')
	output: 
		tumour_vcf=protected('results/variant_analysis/{analysis_type}/{sample}.{nonoverlapping_id}.vcf')
	wildcard_constraints:
		sample='(JBLAB-[0-9]+|IM_[0-9]+)'
	threads: 4
	shell: '../octopus/bin/octopus \
				-C cancer \
				--disable-downsampling \
				--somatics-only \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--min-expected-somatic-frequency 0.03 \
				--min-credible-somatic-frequency 0.01 \
				--max-somatic-haplotypes 2 \
				--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.50 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.90 | SD > 0.90 | FRF > 0.5 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | AD < 1 | AF < 0.0001" \
				--annotations SB SD AF FRF \
				--threads \
				-w temp/ \
				-I {input.tumour_bams} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'

rule filter_octopus_raw_calls_TP53:
	input: 
		filtered_vcf=rules.octopus_tp53.output
	output: 'results/variant_analysis/{analysis_type}/{sample}.{nonoverlapping_id}.filtered.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

# this is less restrictive than expecting matching genotypes
rule remove_homozygous_reference_calls:
	input: rules.filter_octopus_raw_calls_TP53.output
	output: 'results/variant_analysis/{analysis_type}/{sample}.{nonoverlapping_id}.filtered4.vcf'
	shell: 'bash scripts/deprecated/remove_records_with_homozygous_reference_calls.sh {input} > {output}'

rule bgzip_vcf_sample_level:
	input: rules.remove_homozygous_reference_calls.output
	output: temp('results/variant_analysis/{analysis_type}/{sample}.{nonoverlapping_id}.vcf.gz')
	shell: '/home/bioinformatics/software/htslib/htslib-1.6/bin/bgzip < {input} > {output}'

rule index_compressed_vcf_sample_level:
	input: rules.bgzip_vcf_sample_level.output
	output: 'results/variant_analysis/{analysis_type}/{sample}.{nonoverlapping_id}.vcf.gz.csi'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools index {input}'

rule concat_vcfs_sample_level:
	input: 
		compressed_vcfs=lambda wildcards: expand('results/variant_analysis/TP53/{sample}.{nonoverlapping_id}.vcf.gz', nonoverlapping_id=[1,2,3,4,5], sample=wildcards.sample),
		compressed_vcf_indexes=lambda wildcards: expand('results/variant_analysis/TP53/{sample}.{nonoverlapping_id}.vcf.gz.csi', nonoverlapping_id=[1,2,3,4,5], sample=wildcards.sample)
	output: 'results/variant_analysis/TP53/{sample}.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools concat --allow-overlaps {input.compressed_vcfs} -O v -o {output}'
