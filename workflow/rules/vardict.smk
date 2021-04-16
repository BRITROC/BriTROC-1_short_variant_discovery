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

rule convert_picard_interval_to_bed:
	input: 
		targets='resources/panel_28_targets.interval_list',
		amplicons='resources/panel_28_amplicons.interval_list'
	output: 'resources/panel_28.bed'
	script: '../scripts/convert_picard_interval_to_bed.R'

rule vardict:
	input:
		reference_genome=config['reference_genome'],
		path_to_bed_file=rules.convert_picard_interval_to_bed.output,
		tumour_bams=get_tumour_bam_files,
		normal_bams=get_normal_bam_files
	output: 'results/tumour_sample_vcfs_vardict/{sample}.vcf'
	shell: '/scratchb/jblab/bradle02/libraries/VarDict-1.8.2/bin/VarDict \
				-G {input.reference_genome} \
				-f 0.01 \
				-z \
				-U \
				-deldupvar \
				--amplicon 10:0.95 \
				-N tumor_sample_name \
				-b {input.tumour_bams[0]} \
				{input.path_to_bed_file} | \
				/scratchb/jblab/bradle02/libraries/VarDict-1.8.2/bin/teststrandbias.R | \
				/scratchb/jblab/bradle02/libraries/VarDict-1.8.2/bin/var2vcf_valid.pl \
				-N tumour_sample_name \
				-E \
				-f 0.01 > {output}'
