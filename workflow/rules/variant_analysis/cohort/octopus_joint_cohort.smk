def get_tumour_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = all_tumour_metadata[(all_tumour_metadata.fk_britroc_number == int(wildcards.patient_id))]

	bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.place_holder.panel_28.bam', zip, SLX_ID=test_sample_metadata['fk_slx'], barcodes=test_sample_metadata['fk_barcode'], flowcell=test_sample_metadata['flowcell'], lane=test_sample_metadata['lane']) 

	bam_files = []

	for bam_file_name in bam_files_tmp:
		new_bam_name = bam_file_name.replace('place_holder', wildcards.nonoverlapping_id)
		bam_files.append(new_bam_name)

	return(bam_files)

def get_tumour_bam_index_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = all_tumour_metadata[(all_tumour_metadata.fk_britroc_number == int(wildcards.patient_id))]

	bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.place_holder.panel_28.bai', zip, SLX_ID=test_sample_metadata['fk_slx'], barcodes=test_sample_metadata['fk_barcode'], flowcell=test_sample_metadata['flowcell'], lane=test_sample_metadata['lane']) 

	bam_files = []

	for bam_file_name in bam_files_tmp:
		new_bam_name = bam_file_name.replace('place_holder', wildcards.nonoverlapping_id)
		bam_files.append(new_bam_name)

	return(bam_files)

rule convert_bed6_to_oct_format:
	input:  'resources/nonoverlapping.panel_28.targets.{nonoverlapping_id}.bed'                                           
	output: 'resources/panel_28.{nonoverlapping_id}.targets.oct'
	script: '../../../scripts/octopus_formatting/convert_bed6_to_octopus.R'

rule octopus:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.convert_bed6_to_oct_format.output, 
		tumour_bams=get_tumour_bam_files,
		tumour_bam_indexes=get_tumour_bam_index_files
	output: 
		tumour_vcf='results/variant_analysis/cohort/{patient_id}.{nonoverlapping_id}.vcf',
	threads: 4
	wildcard_constraints:
		nonoverlapping_id='[1-9]'
	shell: '../octopus/bin/octopus \
				-C cancer \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--disable-downsampling \
				--forest resources/germline.v0.7.2.forest \
				--somatic-forest resources/somatic.v0.7.2.forest \
				--annotations SB SD AF AD FRF \
				--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.25 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.90 | SD > 0.90 | FRF > 0.5 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | AD < 1 | AF < 0.0001" \
				--min-expected-somatic-frequency 0.03 \
				--min-credible-somatic-frequency 0.01 \
				--threads \
				-w temp/ \
				-I {input.tumour_bams} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'


