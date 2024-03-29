# note: Can't use SelectVariant with octopus output as they are not compatible

# copyright Thomas Bradley 2023 ('thomas.bradley@cruk.cam.ac.uk')

rule filter_octopus_raw_calls:
	input: 
		filtered_vcf=rules.octopus.output
	output: 'results/variant_analysis/matched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered.vcf'
	conda: 'config/bcftools.yaml'
	shell: 'bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

rule filter_calls2:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.convert_bed6_to_oct_format.output,
		tumour_bams=lambda wildcards: get_bam_files(wildcards, 'bam', 'tumour_panel_28'),
		tumour_bam_indexes=lambda wildcards: get_bam_files(wildcards, 'bai', 'tumour_panel_28'),
		normal_bams=lambda wildcards: get_bam_files(wildcards, 'bam', 'normal'),
		normal_bam_indexes = lambda wildcards: get_bam_files(wildcards, 'bai', 'normal'),
		vcf_file=rules.filter_octopus_raw_calls.output
	output: 
		tumour_vcf='results/variant_analysis/matched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered2.vcf',
	threads: 4
	params:
		normal_sample_identifier=get_normal_sample_names
	container: 'docker://dancooke/octopus'
	shell:   'octopus \
				-C cancer \
				--allow-octopus-duplicates \
				--allow-marked-duplicates \
				--disable-downsampling \
				--annotations SB SD AF AD FRF \
				--threads \
				-w temp/ \
				--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.25 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.90 | SD > 0.90 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | FRF > 0.5 | NC > 1 | AD < 1 | AF < 0.03" \
				--filter-vcf {input.vcf_file} \
				-N {params.normal_sample_identifier[0]} \
				-I {input.tumour_bams} {input.normal_bams[0]} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'

rule bgzip_vcf:
	input: rules.filter_calls2.output.tumour_vcf
	output: 
		compressed_vcf=temp('results/variant_analysis/matched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz')
	conda: 'config/htslib.yaml'
	shell: 'bgzip < {input} > {output.compressed_vcf}'

rule index_compressed_vcf:
	input: rules.bgzip_vcf.output
	output: 'results/variant_analysis/matched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz.csi'
	conda: 'config/bcftools.yaml'
	shell: 'bcftools index {input}'

rule concat_vcfs:
	input:
		compressed_vcfs=lambda wildcards: expand('results/variant_analysis/matched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz', analysis_type=wildcards.analysis_type, nonoverlapping_id=get_nonoverlapping_id_list(wildcards), patient_id=wildcards.patient_id),
		compressed_vcf_indexes=lambda wildcards: expand('results/variant_analysis/matched/{analysis_type}/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz.csi', analysis_type=wildcards.analysis_type, nonoverlapping_id=get_nonoverlapping_id_list(wildcards), patient_id=wildcards.patient_id)
	output: 'results/variant_analysis/matched/{analysis_type}/{patient_id}.filtered2.vcf'
	conda: 'config/bcftools.yaml'
	shell: 'bcftools concat --allow-overlaps {input.compressed_vcfs} -O v -o {output}'
