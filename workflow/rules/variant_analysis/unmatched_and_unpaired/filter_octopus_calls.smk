rule filter_unmatched_raw_calls_by_random_forest:
	input: 
		filtered_vcf=rules.octopus_unmatched_with_random_forests.output
	output: 'results/variant_analysis/unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

# Note well the removal of the FRF filter for this particular call to octopus
rule octopus_unmatched_with_hard_filter_annotations:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.convert_bed6_to_oct_format.output,
		tumour_bams=lambda wildcards: get_bam_files(wildcards, 'bam', 'tumour_panel_28'),
		tumour_bam_indexes=lambda wildcards: get_bam_files(wildcards, 'bai', 'tumour_panel_28'),
		vcf_file=rules.filter_unmatched_raw_calls_by_random_forest.output
	output: 
		tumour_vcf='results/variant_analysis/unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered2.vcf',
	threads: 4
	shell:   '../octopus/bin/octopus \
				-C cancer \
				--allow-octopus-duplicates \
				--allow-marked-duplicates \
				--max-somatic-haplotypes 2 \
				--annotations SB SD AF AD FRF \
				--threads \
				-w temp/ \
				--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.25 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.90 | SD > 0.90 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | FRF > 0.5 | AD < 1 | AF < 0.01" \
				--filter-vcf {input.vcf_file} \
				-I {input.tumour_bams} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'

rule bgzip_unmatched_filtered_vcfs:
	input: rules.octopus_unmatched_with_hard_filter_annotations.output.tumour_vcf
	output: 
		compressed_vcf=temp('results/variant_analysis/unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz')
	shell: '/home/bioinformatics/software/htslib/htslib-1.6/bin/bgzip < {input} > {output.compressed_vcf}'

rule index_unmatched_filtered_and_compressed_vcfs:
	input: rules.bgzip_unmatched_filtered_vcfs.output.compressed_vcf
	output: 'results/variant_analysis/unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz.csi'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools index {input}'

rule concat_unmatched_vcfs_across_amplicon_groups:
	input: 
		compressed_vcfs=lambda wildcards: expand('results/variant_analysis/unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz', nonoverlapping_id=get_nonoverlapping_id_list(wildcards), patient_id=wildcards.patient_id, analysis_type=wildcards.analysis_type),
		compressed_vcf_indexes=lambda wildcards: expand('results/variant_analysis/unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz.csi', nonoverlapping_id=get_nonoverlapping_id_list(wildcards), patient_id=wildcards.patient_id, analysis_type=wildcards.analysis_type)
	wildcard_constraints:
		sample='(IM_[0-9]+|JBLAB-[0-9]+)',
		patient_id='[0-9]+'
	output: 'results/variant_analysis/unmatched/{analysis_type}/{patient_id}.filtered2.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools concat --allow-overlaps {input.compressed_vcfs} -O v -o {output}' 
