# note: Can't use SelectVariant with octopus output as they are not compatible

rule filter_octopus_raw_calls:
	input: 
		filtered_vcf=rules.octopus.output
	output: 'results/tumour_sample_vcfs_octopus/{sample}.{nonoverlapping_id}.filtered.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

rule filter_calls2:
	input:
		reference_genome=config['reference_genome'],
		interval_file='resources/panel_6_28.{nonoverlapping_id}.targets.oct',
		tumour_bams=get_tumour_bam_files,
		normal_bams=get_normal_bam_files,
		vcf_file=rules.filter_octopus_raw_calls.output
	output: 
		tumour_vcf='results/tumour_sample_vcfs_octopus/{sample}.{nonoverlapping_id}.filtered2.vcf',
	threads: 4
	params:
		normal_sample_identifier=get_normal_sample_names
	shell:   '../octopus/bin/octopus \
				-C cancer \
				--allow-octopus-duplicates \
				--allow-marked-duplicates \
				--annotations SB SD AF AD FRF \
				--somatics-only \
				--disable-downsampling \
				--threads \
				-w temp/ \
				--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.25 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.90 | SD > 0.90 | FRF > 0.5 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | AD < 1 | AF < 0.03" \
				--filter-vcf {input.vcf_file} \
				-N {params.normal_sample_identifier[0]} \
				-I {input.tumour_bams} {input.normal_bams[0]} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'

rule filter_octopus_hard_filtering:
	input: 
		filtered_vcf='results/tumour_sample_vcfs_octopus/{sample}.filtered2.vcf',
	output: 'results/tumour_sample_vcfs_octopus/{sample}.filtered3.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

rule remove_normal_sample:
	input: 'results/tumour_sample_vcfs_octopus/{sample}.filtered3.vcf'
	output: 'results/tumour_sample_vcfs_octopus/{sample}.filtered4.vcf'
	params: normal_sample_identifier=get_normal_sample_names
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -s "^{params.normal_sample_identifier[0]}" {input} > {output}'

# filters VCF and preserves record in which all samples (except the normal) have the same predicted genotype
rule ensure_matching_genotypes:
	input: 'results/tumour_sample_vcfs_octopus/{sample}.filtered4.vcf'
	output: 'results/tumour_sample_vcfs_octopus/{sample}.filtered5.vcf'
	shell: 'bash workflow/scripts/filter_non_joined_calls/match_genotypes.sh {input} > {output}'

#rule remove_records_with_homozygous_reference_calls:
#	input: 'results/tumour_sample_vcfs_octopus/{sample}.filtered4.vcf'
#	output: 'results/tumour_sample_vcfs_octopus/{sample}.filtered5.vcf'
#	shell: 'bash workflow/scripts/remove_records_with_homozygous_reference_calls.sh {input} > {output}'
