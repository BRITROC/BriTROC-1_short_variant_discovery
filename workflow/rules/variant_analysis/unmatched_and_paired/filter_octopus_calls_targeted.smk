rule bgzip_octopus_unmatched_targeted_vcf:
	input: rules.octopus_unmatched_targeted_calling.output.tumour_vcf
	output: 
		compressed_vcf=temp('results/variant_analysis/unmatched/{analysis_type}/paired/{patient_id}.{nonoverlapping_id}.filtered2.targeted.vcf.gz')
	conda: 'config/htslib.yaml'
	shell: 'bgzip < {input} > {output.compressed_vcf}'

rule index_compressed_octopus_unmatched_vcf_targeted:
	input: rules.bgzip_octopus_unmatched_targeted_vcf.output.compressed_vcf
	output: 'results/variant_analysis/unmatched/{analysis_type}/paired/{patient_id}.{nonoverlapping_id}.filtered2.targeted.vcf.gz.csi'
	conda: 'config/bctools.yaml'
	shell: 'bcftools index {input}'

rule concat_vcfs_octopus_unmatched_targeted:
	input: 
		compressed_vcfs=lambda wildcards: expand('results/variant_analysis/unmatched/{analysis_type}/paired/{patient_id}.{nonoverlapping_id}.filtered2.targeted.vcf.gz', nonoverlapping_id=get_nonoverlapping_id_list(wildcards), patient_id=wildcards.patient_id, analysis_type=wildcards.analysis_type),
		compressed_vcf_indexes=lambda wildcards: expand('results/variant_analysis/unmatched/{analysis_type}/paired/{patient_id}.{nonoverlapping_id}.filtered2.targeted.vcf.gz.csi', nonoverlapping_id=get_nonoverlapping_id_list(wildcards), patient_id=wildcards.patient_id, analysis_type=wildcards.analysis_type)
	wildcard_constraints:
		sample='(IM_[0-9]+|JBLAB-[0-9]+)',
		patient_id='[0-9]+'
	output: 'results/variant_analysis/unmatched/{analysis_type}/paired/{patient_id}.filtered2.targeted.vcf'
	conda: 'config/bcftools.yaml'
	shell: 'bcftools concat --allow-overlaps {input.compressed_vcfs} -O v -o {output}' 
