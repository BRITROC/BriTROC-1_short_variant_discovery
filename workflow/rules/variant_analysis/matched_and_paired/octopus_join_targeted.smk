rule get_interval_file_for_targeted_calling:
	input: rules.ensure_tech_rep_genotypes_match.output.tumour_samples_union
	output: 'results/variant_analysis/matched/{analysis_type}/{patient_id}.interval_file.txt'
	script: '../../../scripts/annotate_variants_joined/get_intervals_for_targeted_calling.R'

rule octopus_targeted_calling:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.get_interval_file_for_targeted_calling.output, 
		tumour_bams=lambda wildcards: get_bam_files(wildcards, 'bam', 'tumour_panel_28'),
		tumour_bam_indexes=lambda wildcards: get_bam_files(wildcards, 'bai', 'tumour_panel_28')
	output: 
		tumour_vcf='results/variant_analysis/matched/{analysis_type}/paired/{patient_id}.{nonoverlapping_id}.targeted.vcf'
	threads: 16
	wildcard_constraints:
		nonoverlapping_id='[1-9]'
	shell: 'octopus \
				-C cancer \
				--disable-downsampling \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--forest resources/germline.v0.7.2.forest \
				--somatic-forest resources/somatic.v0.7.2.forest \
				--somatics-only \
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
