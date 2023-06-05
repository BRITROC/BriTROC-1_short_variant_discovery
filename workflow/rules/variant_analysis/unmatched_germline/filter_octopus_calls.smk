rule filter_octopus_germline_raw_calls:
	input: 
		filtered_vcf=rules.octopus_germline_with_hard_filter_annotation.output
	output: 'results/variant_analysis/germline/octopus_unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

rule filter_octopus_germline_calls2:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.convert_bed6_to_oct_format.output,
		normal_bams= lambda wildcards: get_bam_files(wildcards, 'bam', 'normal'),
		normal_bam_indexes = lambda wildcards: get_bam_files(wildcards, 'bai', 'normal'),
		vcf_file=rules.filter_octopus_germline_raw_calls.output
	output: 'results/variant_analysis/germline/octopus_unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered2.vcf'
	threads: 4
	shell:   '../octopus/bin/octopus \
				--allow-octopus-duplicates \
				--allow-marked-duplicates \
				--disable-downsampling \
				--annotations SB SD AF AD FRF \
				--forest resources/germline.v0.7.2.forest \
				--threads \
				--filter-vcf {input.vcf_file} \
				-I {input.normal_bams} \
				--regions-file {input.interval_file} \
				--output {output} \
				-R {input.reference_genome}'

rule filter_octopus_germline_raw_calls2:
	input: rules.filter_octopus_germline_calls2.output
	output: 'results/variant_analysis/germline/octopus_unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered3.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

rule bgzip_vcf_octopus_germline:
	input: rules.filter_octopus_germline_raw_calls2.output
	output: compressed_vcf=temp('results/variant_analysis/germline/octopus_unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered3.vcf.gz')
	shell: '/home/bioinformatics/software/htslib/htslib-1.6/bin/bgzip < {input} > {output.compressed_vcf}'

rule index_compressed_vcf_octopus_germline:
	input: rules.bgzip_vcf_octopus_germline.output.compressed_vcf
	output: 'results/variant_analysis/germline/octopus_unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered3.vcf.gz.csi'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools index {input}'

rule concat_vcfs_octopus_germline:
	input:
		compressed_vcfs=lambda wildcards: expand('results/variant_analysis/germline/octopus_unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered3.vcf.gz', analysis_type=wildcards.analysis_type, nonoverlapping_id=get_nonoverlapping_id_list(wildcards), patient_id=wildcards.patient_id),
		compressed_vcf_indexes=lambda wildcards: expand('results/variant_analysis/germline/octopus_unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered3.vcf.gz.csi', analysis_type=wildcards.analysis_type, nonoverlapping_id=get_nonoverlapping_id_list(wildcards), patient_id=wildcards.patient_id)
	wildcard_constraints:
		sample='(IM_[0-9]+|JBLAB-[0-9]+)',
		patient_id='[0-9]+'
	output: 'results/variant_analysis/germline/octopus_unmatched/{analysis_type}/merged/{patient_id}_unmatched.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools concat --allow-overlaps {input.compressed_vcfs} -O v -o {output}'

rule ensure_genotypes_match:
	input: rules.concat_vcfs_octopus_germline.output
	output: 'results/variant_analysis/germline/octopus_unmatched/{analysis_type}/merged/{patient_id}.vcf'
	shell: 'workflow/scripts/filter_non_joined_calls/match_genotypes.sh {input} > {output}'
