# octopus only permits the use of one normal sample
rule octopus:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.convert_bed6_to_oct_format.output, 
		tumour_bams= lambda wildcards: get_tumour_bam_files(wildcards, 'bam'),
		tumour_bam_indexes = lambda wildcards: get_tumour_bam_files(wildcards, 'bai'),
		normal_bams= get_normal_bam_files
	output: 
		tumour_vcf=protected('results/variant_analysis/matched/{analysis_type}/{patient_id}.{nonoverlapping_id}.vcf')
	threads: 4
	wildcard_constraints:
		nonoverlapping_id='[1-9]'
	params:
		normal_sample_identifier=get_normal_sample_names
	shell: '../octopus/bin/octopus \
				-C cancer \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--forest resources/germline.v0.7.2.forest \
				--somatic-forest resources/somatic.v0.7.2.forest \
				--disable-downsampling \
				--annotations SB SD AF AD FRF \
				--min-expected-somatic-frequency 0.03 \
				--min-credible-somatic-frequency 0.01 \
				--threads \
				-w temp/ \
				-I {input.tumour_bams} {input.normal_bams[0]} \
				-N {params.normal_sample_identifier[0]} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'
