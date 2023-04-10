def get_tumour_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = matched_and_unpaired_somatic_metadata[(matched_and_unpaired_somatic_metadata.fk_britroc_number == int(wildcards.patient_id))]

	# configure 'analysis_type' string
	if wildcards.analysis_type == 'panel_6_28':
		analysis_type = 'panel_6_28'
	elif wildcards.analysis_type == 'panel_28_only':
		analysis_type = 'antijoined_panel_28_6' 

	bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.nonoverlapping_id_place_holder.analysis_type_place_holder.bam', zip, SLX_ID=test_sample_metadata['fk_slx'], barcodes=test_sample_metadata['fk_barcode'], flowcell=test_sample_metadata['flowcell'], lane=test_sample_metadata['lane'])
	bam_files = []

	for bam_file_name in bam_files_tmp:
		new_bam_name = bam_file_name.replace('nonoverlapping_id_place_holder', wildcards.nonoverlapping_id)
		new_bam_name = new_bam_name.replace('analysis_type_place_holder', analysis_type)
		bam_files.append(new_bam_name)

	return(bam_files)

def get_tumour_bam_index_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = matched_and_unpaired_somatic_metadata[(matched_and_unpaired_somatic_metadata.fk_britroc_number == int(wildcards.patient_id))]

	# configure 'analysis_type' string
	if wildcards.analysis_type == 'panel_6_28':
		analysis_type = 'panel_6_28'
	elif wildcards.analysis_type == 'panel_28_only':
		analysis_type = 'antijoined_panel_28_6' 

	bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.nonoverlapping_id_place_holder.analysis_type_place_holder.bai', zip, SLX_ID=test_sample_metadata['fk_slx'], barcodes=test_sample_metadata['fk_barcode'], flowcell=test_sample_metadata['flowcell'], lane=test_sample_metadata['lane'])
	bam_files = []

	for bam_file_name in bam_files_tmp:
		new_bam_name = bam_file_name.replace('nonoverlapping_id_place_holder', wildcards.nonoverlapping_id)
		new_bam_name = new_bam_name.replace('analysis_type_place_holder', analysis_type)
		bam_files.append(new_bam_name)

	return(bam_files)

def get_normal_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = matched_and_unpaired_somatic_metadata[(matched_and_unpaired_somatic_metadata.fk_britroc_number == int(wildcards.patient_id))]
	britroc_number = test_sample_metadata.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	normal_metadata = matched_and_unpaired_germline_metadata[matched_and_unpaired_germline_metadata['fk_britroc_number'] == britroc_number[0]]

	# configure 'analysis_type' string
	if wildcards.analysis_type == 'panel_6_28':
		analysis_type = 'panel_6_28'
	elif wildcards.analysis_type == 'panel_28_only':
		analysis_type = 'antijoined_panel_28_6' 

	bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.nonoverlapping_id_place_holder.analysis_type_place_holder.bam', zip, SLX_ID=normal_metadata['fk_slx'], barcodes=normal_metadata['fk_barcode'], flowcell=normal_metadata['flowcell'], lane=normal_metadata['lane']) 
	bam_files = []

	for bam_file_name in bam_files_tmp:
		new_bam_name = bam_file_name.replace('nonoverlapping_id_place_holder', wildcards.nonoverlapping_id)
		new_bam_name = new_bam_name.replace('analysis_type_place_holder', analysis_type)
		bam_files.append(new_bam_name)

	return(bam_files)

def get_normal_sample_names(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = matched_and_unpaired_somatic_metadata[(matched_and_unpaired_somatic_metadata.fk_britroc_number == int(wildcards.patient_id))]
	britroc_number = test_sample_metadata.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	normal_metadata = matched_and_unpaired_germline_metadata[matched_and_unpaired_germline_metadata['fk_britroc_number'] == britroc_number[0]]

	SLX = normal_metadata.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	
	samples = normal_metadata.set_index('fk_sample', drop=False)
	samples = samples.index.unique().tolist()	

	if SLX[0] in ['SLX-14363']:
		samples.append(samples[0] + '_d')
		samples[0] = '{}.{}'.format(SLX[0], barcodes[0])
		samples[1] = '{}.{}'.format(SLX[0], barcodes[1])
	elif SLX[0] == 'SLX-16247' and samples[0] in ['IM_429','IM_432','IM_437','IM_438']:
		samples.append(samples[0])
		samples[1] = '{}_{}'.format(samples[0], 'd2')
		samples[0] = '{}_{}'.format(samples[0], 'd1')
	#elif SLX[0] == 'SLX-9856':
	#	samples = [sample.replace('-', '') for sample in samples]
	else:
		samples.append(samples[0] + '_d')

	print(samples)

	return(samples)

def get_relevant_bed_file(wildcards):
	if wildcards.analysis_type == 'panel_6_28':
		output_file='resources/panel_6_28.nonoverlapping.targets.{}.bed'.format(wildcards.nonoverlapping_id)
	elif wildcards.analysis_type == 'panel_28_only':
		output_file='resources/nonoverlapping.antijoined_panel_28_6.targets.{}.bed'.format(wildcards.nonoverlapping_id)
	
	return(output_file)

rule convert_bed6_to_oct_format:
	input:  get_relevant_bed_file                                           
	output: 'resources/{analysis_type}.{nonoverlapping_id}.targets.oct'
	script: '../../../scripts/octopus_formatting/convert_bed6_to_octopus.R'

# octopus only permits the use of one normal sample
rule octopus:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.convert_bed6_to_oct_format.output, 
		tumour_bams=get_tumour_bam_files,
		normal_bams=get_normal_bam_files
	output: 
		tumour_vcf=protected('results/variant_analysis/matched/{analysis_type}/{patient_id}.{nonoverlapping_id}.vcf')
	threads: 4
	wildcard_constraints:
		nonoverlapping_id='[1-9]'
	params:
		normal_sample_identifier=get_normal_sample_names
	shell: '../octopus/bin/octopus \
				-C cancer \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--forest resources/germline.v0.7.2.forest \
				--somatic-forest resources/somatic.v0.7.2.forest \
				--disable-downsampling \
				--annotations SB SD AF AD FRF \
				--min-expected-somatic-frequency 0.03 \
				--min-credible-somatic-frequency 0.01 \
				--threads \
				-w temp/ \
				-I {input.tumour_bams} {input.normal_bams[0]} \
				-N {params.normal_sample_identifier[0]} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'
