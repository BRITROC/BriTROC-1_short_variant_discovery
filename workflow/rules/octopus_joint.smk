def get_tumour_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = matched_somatic_metadata[(matched_somatic_metadata.fk_britroc_number == int(wildcards.patient_id))]

	bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.place_holder.panel_6_28.bam', zip, SLX_ID=test_sample_metadata['fk_slx'], barcodes=test_sample_metadata['fk_barcode'], flowcell=test_sample_metadata['flowcell'], lane=test_sample_metadata['lane']) 

	bam_files = []

	for bam_file_name in bam_files_tmp:
		new_bam_name = bam_file_name.replace('place_holder', wildcards.nonoverlapping_id)
		bam_files.append(new_bam_name)

	return(bam_files)

def get_normal_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = matched_somatic_metadata[(matched_somatic_metadata.fk_britroc_number == int(wildcards.patient_id))]
	britroc_number = test_sample_metadata.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	normal_metadata = matched_germline_metadata[matched_germline_metadata['fk_britroc_number'] == britroc_number[0]]

	bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.place_holder.panel_6_28.bam', zip, SLX_ID=normal_metadata['fk_slx'], barcodes=normal_metadata['fk_barcode'], flowcell=normal_metadata['flowcell'], lane=normal_metadata['lane']) 

	bam_files = []

	for bam_file_name in bam_files_tmp:
		new_bam_name = bam_file_name.replace('place_holder', wildcards.nonoverlapping_id)
		bam_files.append(new_bam_name)

	return(bam_files)

def get_normal_sample_names(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = matched_somatic_metadata[(matched_somatic_metadata.fk_britroc_number == int(wildcards.patient_id))]
	britroc_number = test_sample_metadata.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	normal_metadata = matched_germline_metadata[matched_germline_metadata['fk_britroc_number'] == britroc_number[0]]

	SLX = normal_metadata.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	
	samples = normal_metadata.set_index('fk_sample', drop=False)
	samples = samples.index.unique().tolist()	

	if SLX[0] == 'SLX-9629':
		samples[0] = samples[0].replace('-','')
	elif SLX[0] in ['SLX-13630','SLX-11347','SLX-11111','SLX-11110','SLX-11109','SLX-9856']:
		# this is accounting for a bug in the formatting of BAM headers
		samples.append(samples[0] + '_d')
		samples[0] = samples[0].replace('-','')
		samples[1] = samples[1].replace('-','')
	elif SLX[0] == 'SLX-14363':
		samples.append(samples[0] + '_d')
		samples[0] = '{}.{}'.format(SLX[0], barcodes[0])
		samples[1] = '{}.{}'.format(SLX[0], barcodes[1])
	elif SLX[0] == 'SLX-16247' and samples[0] in ['IM_429','IM_432','IM_437','IM_438']:
		samples.append(samples[0])
		samples[1] = '{}_{}'.format(samples[0], 'd2')
		samples[0] = '{}_{}'.format(samples[0], 'd1')
	else:
		samples.append(samples[0] + '_d')

	return(samples)

rule convert_bed6_to_oct_format:
	input:  'resources/panel_6_28.nonoverlapping.targets.{nonoverlapping_id}.bed'                                           
	output: 'resources/panel_6_28.{nonoverlapping_id}.targets.oct'
	script: '../scripts/convert_bed6_to_octopus.R'

rule octopus:
	input:
		reference_genome=config['reference_genome'],
		interval_file='resources/panel_6_28.{nonoverlapping_id}.targets.oct', 
		tumour_bams=get_tumour_bam_files,
		normal_bams=get_normal_bam_files
	output: 
		tumour_vcf='results/tumour_sample_vcfs_octopus/{patient_id}.{nonoverlapping_id}.vcf',
	threads: 4
	wildcard_constraints:
		nonoverlapping_id='[1-9]'
	params:
		normal_sample_identifier=get_normal_sample_names
	shell: '../octopus/bin/octopus \
				-C cancer \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--somatic-forest-model resources/somatic.v0.7.2.forest \
				--somatics-only \
				--disable-downsampling \
				--annotations SB SD AF AD FRF \
				--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.25 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.90 | SD > 0.90 | FRF > 0.5 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | AD < 1 | AF < 0.0001" \
				--min-expected-somatic-frequency 0.03 \
				--min-credible-somatic-frequency 0.01 \
				--threads \
				-w temp/ \
				-I {input.tumour_bams} {input.normal_bams[0]} \
				-N {params.normal_sample_identifier[0]} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'

rule bgzip_vcf:
	input: 'results/tumour_sample_vcfs_octopus/{patient_id}.{nonoverlapping_id}.filtered2.vcf'
	output: 
		compressed_vcf=temp('results/tumour_sample_vcfs_octopus/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz')
	shell: '/home/bioinformatics/software/htslib/htslib-1.6/bin/bgzip < {input} > {output.compressed_vcf}'

rule index_compressed_vcf:
	input: 'results/tumour_sample_vcfs_octopus/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz'
	output: 'results/tumour_sample_vcfs_octopus/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz.csi'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools index {input}'

rule concat_vcfs:
	input: 
		compressed_vcfs=lambda wildcards: expand('results/tumour_sample_vcfs_octopus/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz', nonoverlapping_id=[1,2,3,4], patient_id=wildcards.patient_id),
		compressed_vcf_indexes=lambda wildcards: expand('results/tumour_sample_vcfs_octopus/{patient_id}.{nonoverlapping_id}.filtered2.vcf.gz.csi', nonoverlapping_id=[1,2,3,4], patient_id=wildcards.patient_id)
	wildcard_constraints:
		sample='(IM_[0-9]+|JBLAB-[0-9]+)',
		patient_id='[0-9]+'
	output: 'results/tumour_sample_vcfs_octopus/{patient_id}.filtered2.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools concat --allow-overlaps {input.compressed_vcfs} -O v -o {output}' 


