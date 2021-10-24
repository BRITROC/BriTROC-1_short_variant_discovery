# note: Can't use SelectVariant with octopus output as they are not compatible

rule filter_octopus_raw_calls_TP53:
	input: 
		filtered_vcf=rules.octopus_tp53.output
	output: 'results/variant_analysis/TP53/{sample}.{nonoverlapping_id}.filtered.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

rule filter_calls2_TP53:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.convert_bed6_to_oct_format_tp53.output,
		tumour_bams=get_tumour_bam_files,
		vcf_file=rules.filter_octopus_raw_calls_TP53.output
	output: 
		tumour_vcf='results/variant_analysis/TP53/{sample}.{nonoverlapping_id}.filtered2.vcf',
	threads: 4
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
				-I {input.tumour_bams} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'

rule filter_octopus_hard_filtering:
	input: filtered_vcf=rules.filter_calls2_TP53.output,
	output: 'results/variant_analysis/TP53/{sample}.filtered3.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

# filters VCF and preserves record in which all samples (except the normal) have the same predicted genotype
rule ensure_matching_genotypes:
	input: rules.filter_octopus_hard_filtering.output
	output: 'results/variant_analysis/TP53/{sample}.filtered4.vcf'
	shell: 'bash workflow/scripts/filter_non_joined_calls/match_genotypes.sh {input} > {output}'
