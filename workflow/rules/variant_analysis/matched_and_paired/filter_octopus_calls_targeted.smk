rule filter_octopus_raw_calls_targeted:
	input: 
		filtered_vcf=rules.octopus_targeted_calling.output
	output: 'results/variant_analysis/matched/{analysis_type}/paired/{patient_id}.{nonoverlapping_id}.filtered.targeted.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

rule bgzip_vcf_targeted:
	input: rules.filter_octopus_raw_calls_targeted.output
	output: 
		compressed_vcf=temp('results/variant_analysis/matched/{analysis_type}/paired/{patient_id}.{nonoverlapping_id}.filtered2.targeted.vcf.gz')
	shell: '/home/bioinformatics/software/htslib/htslib-1.6/bin/bgzip < {input} > {output.compressed_vcf}'

rule index_compressed_vcf_targeted:
	input: rules.bgzip_vcf_targeted.output.compressed_vcf
	output: 'results/variant_analysis/matched/{analysis_type}/paired/{patient_id}.{nonoverlapping_id}.filtered2.targeted.vcf.gz.csi'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools index {input}'

rule concat_vcfs_targeted:
	input: 
		compressed_vcfs=lambda wildcards: expand('results/variant_analysis/matched/{analysis_type}/paired/{patient_id}.{nonoverlapping_id}.filtered2.targeted.vcf.gz', get_nonoverlapping_id_list(wildcards), patient_id=wildcards.patient_id, analysis_type=wildcards.analysis_type),
		compressed_vcf_indexes=lambda wildcards: expand('results/variant_analysis/matched/{analysis_type}/paired/{patient_id}.{nonoverlapping_id}.filtered2.targeted.vcf.gz.csi', get_nonoverlapping_id_list(wildcards), patient_id=wildcards.patient_id, analysis_type=wildcards.analysis_type)
	output: 'results/variant_analysis/matched/{analysis_type}/paired/{patient_id}.filtered2.targeted.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools concat --allow-overlaps {input.compressed_vcfs} -O v -o {output}'
