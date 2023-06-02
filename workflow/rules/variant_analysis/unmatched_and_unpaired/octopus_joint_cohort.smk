# NB: we use both germline and somatic random forests as we are calling both variant types even though there is no matched normal
rule octopus_unmatched_with_random_forests:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.convert_bed6_to_oct_format.output, 
		tumour_bams= lambda wildcards: get_bam_files(wildcards, 'bam', 'tumour_panel_28'),
		tumour_bam_indexes= lambda wildcards: get_bam_files(wildcards, 'bai', 'tumour_panel_28')
	output: 
		tumour_vcf=protected('results/variant_analysis/{matched_or_unmatched}/{analysis_type}/{patient_id}.{nonoverlapping_id}.vcf')
	threads: 16
	wildcard_constraints:
		nonoverlapping_id='[1-9]'
	shell: '../octopus/bin/octopus \
				-C cancer \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--disable-downsampling \
				--forest resources/germline.v0.7.2.forest \
				--somatic-forest resources/somatic.v0.7.2.forest \
				--max-somatic-haplotypes 2 \
				--annotations SB SD AF AD FRF \
				--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.25 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.90 | SD > 0.90 | FRF > 0.5 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | AD < 1 | AF < 0.0001" \
				--min-expected-somatic-frequency 0.03 \
				--min-credible-somatic-frequency 0.01 \
				--threads \
				-w temp/ \
				-I {input.tumour_bams} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'
