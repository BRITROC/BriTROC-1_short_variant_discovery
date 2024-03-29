rule get_interval_file_for_targeted_calling_for_unmatched_analysis:
	input: rules.ensure_tech_rep_genotypes_match_for_unmatched_variants.output.tumour_samples_union
	output: 'results/variant_analysis/unmatched/{analysis_type}/paired/{patient_id}.interval_file.txt'
	script: '../../../scripts/annotate_variants_joined/get_intervals_for_targeted_calling.R'

# rerun octopus though this time with more liberal hard thresholds
rule octopus_unmatched_targeted_calling:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.get_interval_file_for_targeted_calling_for_unmatched_analysis.output, 
		tumour_bams=lambda wildcards: get_bam_files(wildcards, 'bam', 'tumour_panel_28'),
		tumour_bam_indexes=lambda wildcards: get_bam_files(wildcards, 'bai', 'tumour_panel_28'),
	output: 
		tumour_vcf='results/variant_analysis/unmatched/{analysis_type}/paired/{patient_id}.{nonoverlapping_id}.targeted.vcf',
	threads: 16
	wildcard_constraints:
		nonoverlapping_id='[1-9]'
	container: 'docker://dancooke/octopus'
	shell: 'octopus \
				-C cancer \
				--disable-downsampling \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--max-somatic-haplotypes 2 \
				--annotations SB SD AF AD FRF \
				--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.001 | AFB > 0.50 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.90 | SD > 0.90 | FRF > 0.5 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | AD < 1 | AF < 0.0001" \
				--min-expected-somatic-frequency 0.03 \
				--min-credible-somatic-frequency 0.01 \
				--threads \
				-w temp/ \
				-I {input.tumour_bams} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'
