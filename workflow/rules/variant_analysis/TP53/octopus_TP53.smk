def get_tumour_bam_files(wildcards):
	test_sample_metadata = somatic_tp53_metadata[(somatic_tp53_metadata.fk_sample == wildcards.sample)]

	bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.place_holder.union_of_tp53.bam', zip, SLX_ID=test_sample_metadata['fk_slx'], barcodes=test_sample_metadata['fk_barcode'], flowcell=test_sample_metadata['flowcell'], lane=test_sample_metadata['lane']) 

	bam_files = []

	for bam_file_name in bam_files_tmp:
		new_bam_name = bam_file_name.replace('place_holder', wildcards.nonoverlapping_id)
		bam_files.append(new_bam_name)

	return(bam_files)

rule convert_bed6_to_oct_format_tp53:
	input:  rules.create_nonoverlapping_amplicons_TP53.output.nonoverlapping_target_intervals                                     
	output: 'resources/union_of_tp53.{nonoverlapping_id}.targets.oct'
	script: '../scripts/octopus_formatting/convert_bed6_to_octopus.R'

rule octopus_tp53:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.convert_bed6_to_oct_format_tp53.output,
		tumour_bams=get_tumour_bam_files,
	output: 
		tumour_vcf='results/variant_analysis/TP53/{sample}.{nonoverlapping_id}.vcf'
	wildcard_constraints:
		sample='(JBLAB-[0-9]+|IM_[0-9]+)'
	threads: 4
	shell: '../octopus/bin/octopus \
				-C cancer \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--min-expected-somatic-frequency 0.03 \
				--min-credible-somatic-frequency 0.01 \
				--max-somatic-haplotypes 2 \
				--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.0001 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.90 | SD > 0.90 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | AD < 1 | AF < 0.0001" \
				--annotations SB SD AF FRF \
				--threads \
				-w temp/ \
				-I {input.tumour_bams} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'
rule bgzip_vcf_sample_level:
	input: rules.octopus_tp53.output
	output: temp('results/variant_analysis/TP53/{sample}.{nonoverlapping_id}.vcf.gz')
	wildcard_constraints:
		sample='(JBLAB-[0-9]+|IM_[0-9]+)'
	shell: '/home/bioinformatics/software/htslib/htslib-1.6/bin/bgzip < {input} > {output.compressed_vcf}'

rule index_compressed_vcf_sample_level:
	input: rules.bgzip_vcf_sample_level.output
	output: 'results/variant_analysis/TP53/{sample}.{nonoverlapping_id}.vcf.gz.csi'
	wildcard_constraints:
		sample='(JBLAB-[0-9]+|IM_[0-9]+)'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools index {input}'

rule concat_vcfs_sample_level:
	input: 
		compressed_vcfs=lambda wildcards: expand('results/variant_analysis/{sample}.{nonoverlapping_id}.vcf.gz', nonoverlapping_id=[1,2,3,4,5], sample=wildcards.sample),
		compressed_vcf_indexes=lambda wildcards: expand('results/variant_analysis/{sample}.{nonoverlapping_id}.vcf.gz.csi', nonoverlapping_id=[1,2,3,4,5], sample=wildcards.sample)
	output: 'results/variant_analysis/TP53/{sample}.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools concat --allow-overlaps {input.compressed_vcfs} -O v -o {output}'