# NB: we use both germline and somatic random forests as we are calling both variant types even though there is no matched normal
rule octopus_unmatched_with_random_forests:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.convert_bed6_to_oct_format.output, 
		tumour_bams= lambda wildcards: get_bam_files(wildcards, 'bam', 'tumour_panel_28'),
		tumour_bam_indexes= lambda wildcards: get_bam_files(wildcards, 'bai', 'tumour_panel_28')
	output: 
		tumour_vcf=protected('results/variant_analysis/unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.vcf')
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
				--min-expected-somatic-frequency 0.03 \
				--min-credible-somatic-frequency 0.01 \
				--threads \
				-w temp/ \
				-I {input.tumour_bams} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'
