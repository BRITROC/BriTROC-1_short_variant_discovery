def get_tumour_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = somatic_tp53_metadata[(somatic_tp53_metadata.fk_sample == wildcards.sample)]

	#SLX = test_sample_metadata.set_index('fk_slx', drop=False)
	#SLX = SLX.index.unique().tolist()
	#barcodes = test_sample_metadata.set_index('fk_barcode', drop=False)
	#barcodes = barcodes.index.unique().tolist()
	#lane = test_sample_metadata.set_index('lane', drop=False)
	#lane = lane.index.unique().tolist()
	#flowcell = test_sample_metadata.set_index('fk_run', drop=False)
	#flowcell = flowcell.index.unique().tolist()
	#flowcell = flowcell[0].split('_')

	bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.place_holder.union_of_tp53.bam', zip, SLX_ID=test_sample_metadata['fk_slx'], barcodes=test_sample_metadata['fk_barcode'], flowcell=test_sample_metadata['flowcell'], lane=test_sample_metadata['lane']) 

	bam_files = []

	for bam_file_name in bam_files_tmp:
		new_bam_name = bam_file_name.replace('place_holder', wildcards.nonoverlapping_id)
		bam_files.append(new_bam_name)

	return(bam_files)

#rule convert_interval_list_to_bed_format:
#	input: 
#		targets='resources/{amplicon_panel_id}_amplicons.targets.interval_list',
#		amplicons='resources/{amplicon_panel_id}_amplicons.amplicons.interval_list'
#	output: 'resources/{amplicon_panel_id}.bed'
#	script: '../scripts/convert_picard_interval_to_bed.R'

#rule convert_bed_to_oct_format:
#	input:  'tp53.nonoverlapping.targets.{nonoverlapping_id}.bed'                                           #rules.convert_interval_list_to_bed_format.output
#	output: 'resources/union_of_tp53.{nonoverlapping_id}.targets.oct'
#	script: '../scripts/convert_bed_to_octopus.R'

rule convert_bed6_to_oct_format:
	input:  'tp53.nonoverlapping.targets.{nonoverlapping_id}.bed'                                           #rules.convert_interval_list_to_bed_format.output
	output: 'resources/union_of_tp53.{nonoverlapping_id}.targets.oct'
	script: '../scripts/convert_bed6_to_octopus.R'

rule octopus:
	input:
		reference_genome=config['reference_genome'],
		interval_file='resources/union_of_tp53.{nonoverlapping_id}.targets.oct',             #'resources/intersected_panel_6_28.oct',
		tumour_bams=get_tumour_bam_files,
	output: 
		tumour_vcf='results/tumour_sample_vcfs_octopus/{sample}.{nonoverlapping_id}.vcf'
		#realigned_bam=directory('results/realigned_bams/{sample}_octopus')
	wildcard_constraints:
		sample='(JBLAB-[0-9]+|IM_[0-9]+)'
	threads: 4
	shell: '../octopus/bin/octopus \
				-C cancer \
				--disable-downsampling \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--min-expected-somatic-frequency 0.05 \
				--min-credible-somatic-frequency 0.01 \
				--max-somatic-haplotypes 2 \
				--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.0001 | AFB > 0.95 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.90 | SD > 0.90 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | AD < 1 | AF < 0.0001" \
				--annotations SB SD AF FRF \
				--threads \
				-w temp/ \
				-I {input.tumour_bams} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'
rule bgzip_vcf:
	input: 'results/tumour_sample_vcfs_octopus/{sample}.{nonoverlapping_id}.vcf'
	output: 
		compressed_vcf=temp('results/tumour_sample_vcfs_octopus/{sample}.{nonoverlapping_id}.vcf.gz')
	shell: '/home/bioinformatics/software/htslib/htslib-1.6/bin/bgzip < {input} > {output.compressed_vcf}'

rule index_compressed_vcf:
	input: 'results/tumour_sample_vcfs_octopus/{sample}.{nonoverlapping_id}.vcf.gz'
	output: 'results/tumour_sample_vcfs_octopus/{sample}.{nonoverlapping_id}.vcf.gz.csi'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools index {input}'

rule concat_vcfs:
	input: 
		compressed_vcfs=lambda wildcards: expand('results/tumour_sample_vcfs_octopus/{sample}.{nonoverlapping_id}.vcf.gz', nonoverlapping_id=[1,2,3,4,5], sample=wildcards.sample),
		compressed_vcf_indexes=lambda wildcards: expand('results/tumour_sample_vcfs_octopus/{sample}.{nonoverlapping_id}.vcf.gz.csi', nonoverlapping_id=[1,2,3,4,5], sample=wildcards.sample)
	output: 'results/tumour_sample_vcfs_octopus/{sample}.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools concat --allow-overlaps {input.compressed_vcfs} -O v -o {output}' 

#rule merge_vcfs:
#	input: lambda wildcards: expand('results/tumour_sample_vcfs_octopus/{sample}.{nonoverlapping_id}.vcf', nonoverlapping_id=[1,2,3,4,5], sample=wildcards.sample),
#	output: 'results/tumour_sample_vcfs_octopus/{sample}.vcf'
#	shell: 'cat {input[0]} > output_tmp & \
		#			grep -v "#" {input[1]} >> output_tmp & \
		#	grep -v "#" {input[2]} >> output_tmp & \
		#	grep -v "#" {input[3]} >> output_tmp & \
		#	grep -v "#" {input[4]} >> output_tmp & \
		#	mv output_tmp {output}'

#rule filter_octopus_calls:
#	input:
#		reference_genome=config['reference_genome'],
#		interval_file='resources/union_of_tp53.oct',                                #'resources/intersected_panel_6_28.oct',
#		tumour_bams=get_tumour_bam_files,
#		vcf_file='results/tumour_sample_vcfs_octopus/{sample}.filtered.vcf'
#	output: 
#		tumour_vcf='results/tumour_sample_vcfs_octopus/{sample}.filtered2.vcf',
#	threads: 4
#	run:
#			shell ('../octopus/bin/octopus \
		#						-C cancer \
		#		--allow-octopus-duplicates \
		#		--threads \
		#		-w temp/ \
		#		--filter-vcf {input.vcf_file} \
		#		--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 1.20 | SD > 1.20 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | AD < 1 | AF < 0.0001" \
		#		--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.0001 | AFB > 0.95 | SB > 1.20 | BQ < 15 | DP < 1 | ADP < 1" \
		#		-I {input.tumour_bams} \
		#		--regions-file {input.interval_file} \
		#		--output {output.tumour_vcf} \
		#		-R {input.reference_genome}')


#rule filter_octopus_calls2:
#	input: 'results/tumour_sample_vcfs_octopus/{sample}.filtered3.vcf'
#	output: 'results/tumour_sample_vcfs_octopus/{sample}.filtered4.vcf'
#	shell: "grep -E \(##\|HSS\|#\) {input} > {output}"

#rule filter_octopus_calls2:
#	input:
#		reference_genome=config['reference_genome'],
#		interval_file='resources/union_of_tp53.oct',                             #'resources/intersected_panel_6_28.oct',
#		tumour_bams=get_tumour_bam_files,
#		vcf_file='results/tumour_sample_vcfs_octopus/{sample}.filtered.vcf'
#	output: 
#		tumour_vcf='results/tumour_sample_vcfs_octopus/{sample}.filtered2.vcf',
#	threads: 4
#	run:
#			shell ('../octopus/bin/octopus \
		#				-C cancer \
		#				--allow-octopus-duplicates \
		#		--threads \
		#		-w temp/ \
		#		--filter-vcf {input.vcf_file} \
		#		--somatic-forest resources/somatic.v0.7.2.forest \
		#		--somatics-only \
		#		-I {input.tumour_bams} \
		#		--regions-file {input.interval_file} \
		#		--output {output.tumour_vcf} \
		#		-R {input.reference_genome}')

#rule analyse_tp53:
#	input: rules.octopus.output.tumour_vcf
#	output: 'results/tumour_sample_vcfs_octopus/{sample}_tp53.txt'
#	script: '../scripts/process_oct_output.R'
