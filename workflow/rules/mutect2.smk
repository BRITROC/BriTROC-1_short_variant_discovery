#!/bin/env python

#rule all:
#	input: 
		#expand('tumour_vcfs/{SLX}_{barcode}_{flowcell}_{lane}_pileups.table', zip, SLX=metadata['fk_slx'], barcode=metadata['fk_barcode'], flowcell=metadata['flowcell'], lane=metadata['lane']),
		#expand('tumour_vcfs/{SLX}_{barcode}_{flowcell}_{lane}_contamination.table', zip, SLX=metadata['fk_slx'], barcode=metadata['fk_barcode'], flowcell=metadata['flowcell'], lane=metadata['lane']),
		#expand('tumour_vcfs/{SLX}_{barcode}_{flowcell}_{lane}_segmentation.table', zip, SLX=metadata['fk_slx'], barcode=metadata['fk_barcode'], flowcell=metadata['flowcell'], lane=metadata['lane']),

def get_tumour_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = paired_somatic_metadata[(paired_somatic_metadata.fk_sample == wildcards.sample)].head(n=2)

	SLX = test_sample_metadata.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	barcodes = test_sample_metadata.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()
	lane = test_sample_metadata.set_index('lane', drop=False)
	lane = lane.index.unique().tolist()
	flowcell = test_sample_metadata.set_index('fk_run', drop=False)
	flowcell = flowcell.index.unique().tolist()
	flowcell = flowcell[0].split('_')

	bam_file_1 = '../SLX/{}/bam/{}.{}.{}.s_{}.bam'.format(SLX[0],  ''.join(SLX), barcodes[0], flowcell[-1], lane[0])
	bam_file_2 = '../SLX/{}/bam/{}.{}.{}.s_{}.bam'.format(SLX[0],  ''.join(SLX), barcodes[1], flowcell[-1], lane[0])

	return([bam_file_1, bam_file_2])

def get_normal_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = paired_somatic_metadata[(paired_somatic_metadata.fk_sample == wildcards.sample)].head(n=2)
	britroc_number = test_sample_metadata.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	patient_metadata = paired_somatic_metadata[paired_somatic_metadata['fk_britroc_number'] == britroc_number[0]]

	test_tumor_type = test_sample_metadata.set_index('type', drop=False)
	test_tumour_type = test_tumor_type.index.unique().tolist()
	
	if test_tumour_type[0] == 'relapse': 
		reference_tumour_metadata = patient_metadata[patient_metadata['type'] == 'archival'].head(n=2)
	elif test_tumour_type[0] == 'archival':
		reference_tumour_metadata = patient_metadata[patient_metadata['type'] == 'relapse'].head(n=2)

	SLX = reference_tumour_metadata.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()

	barcodes = reference_tumour_metadata.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()

	lane = reference_tumour_metadata.set_index('lane', drop=False)
	lane = lane.index.unique().tolist()
	
	flowcell = reference_tumour_metadata.set_index('fk_run', drop=False)
	flowcell = flowcell.index.unique().tolist()
	flowcell = flowcell[0].split('_')

	if SLX[0] == 'SLX-9629':

		bam_file_1 = '../SLX/{}/samplebam/{}.{}.bam'.format(SLX[0], SLX[0], barcodes[0])
		return([bam_file_1])

	else:
		bam_file_1 = '../SLX/{}/bam/{}.{}.{}.s_{}.bam'.format(SLX[0],  ''.join(SLX), barcodes[0], flowcell[-1], lane[0])
		bam_file_2 = '../SLX/{}/bam/{}.{}.{}.s_{}.bam'.format(SLX[0],  ''.join(SLX), barcodes[1], flowcell[-1], lane[0])

		return([bam_file_1, bam_file_2])

def get_normal_sample_names(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = paired_somatic_metadata[(paired_somatic_metadata.fk_sample == wildcards.sample)].head(n=2)
	britroc_number = test_sample_metadata.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	patient_metadata = paired_somatic_metadata[paired_somatic_metadata['fk_britroc_number'] == britroc_number[0]]

	test_tumor_type = test_sample_metadata.set_index('type', drop=False)
	test_tumour_type = test_tumor_type.index.unique().tolist()
	
	if test_tumour_type[0] == 'relapse': 
		reference_tumour_metadata = patient_metadata[patient_metadata['type'] == 'archival'].head(n=2)
	elif test_tumour_type[0] == 'archival':
		reference_tumour_metadata = patient_metadata[patient_metadata['type'] == 'relapse'].head(n=2)

	SLX = reference_tumour_metadata.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	
	samples = reference_tumour_metadata.set_index('fk_sample', drop=False)
	samples = samples.index.unique().tolist()	

	barcodes = reference_tumour_metadata.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()

	if SLX[0] == 'SLX-9629':
		pass
		#samples[0] = samples[0].replace('-','')
	elif SLX[0] == 'SLX-13630':
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

rule mutect2:
	input:
		reference_genome=config['reference_genome'],
		panel_of_normals=rules.create_panel_of_normals.output,
		germline_resource='resources/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf.bgz',
		interval_file='resources/panel_28.interval_list',
		tumour_bams=get_tumour_bam_files,
		normal_bams=get_normal_bam_files
	output: 
		tumour_vcf='results/tumour_sample_vcfs/{sample}.vcf',
		f1r2='results/tumour_sample_vcfs/{sample}_f1r2.tar.gz',
		bam_output='results/tumour_sample_vcfs/{sample}.bam'
	threads: 4
	params:
		normal_sample_identifier=get_normal_sample_names
	run:
		print(input.normal_bams)
		if len(params.normal_sample_identifier) == 2:
			print('2 normal sample bam files')
			shell ('/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk Mutect2 \
				--input {input.tumour_bams[0]} \
				--input {input.tumour_bams[1]} \
				--input {input.normal_bams[0]} \
				--input {input.normal_bams[1]} \
				--normal-sample {params.normal_sample_identifier[0]} \
				--normal-sample {params.normal_sample_identifier[1]} \
				--intervals {input.interval_file} \
				--panel-of-normals {input.panel_of_normals} \
				--germline-resource {input.germline_resource} \
				--output {output.tumour_vcf} \
				--f1r2-tar-gz {output.f1r2} \
				--minimum-allele-fraction 0.00 \
				--bam-output {output.bam_output} \
				--reference {input.reference_genome}')
		else:
			print('1 normal sample bam file')
			shell ('/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk Mutect2 \
				--input {input.tumour_bams[0]} \
				--input {input.tumour_bams[1]} \
				--input {input.normal_bams[0]} \
				--normal-sample {params.normal_sample_identifier[0]} \
				--intervals {input.interval_file} \
				--panel-of-normals {input.panel_of_normals} \
				--germline-resource {input.germline_resource} \
				--output {output.tumour_vcf} \
				--f1r2-tar-gz {output.f1r2} \
				--minimum-allele-fraction 0.00 \
				--bam-output {output.bam_output} \
				--reference {input.reference_genome}')
