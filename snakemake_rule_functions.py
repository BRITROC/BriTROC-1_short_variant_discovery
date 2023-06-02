# a file containing functions used within snakemake rules
# This is distinct from the list of functions used within R scripts

def get_bam_files(wildcards, bam_or_bai, sample_type):

	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	# NB: the paired metadata table is a subset of the unpaired table - therefore using the unpaired table for the paired analysis is not inappropriate

	if sample_type == 'tumour_panel_28':
		test_sample_metadata = panel_28_tumour_sequencing_metadata[(panel_28_tumour_sequencing_metadata.fk_britroc_number == int(wildcards.patient_id))]
	if sample_type == 'tumour_TP53':
		test_sample_metadata = TP53_tumour_sequencing_metadata[(TP53_tumour_sequencing_metadata.fk_sample == wildcards.sample)]
	elif sample_type == 'normal':
		test_sample_metadata = nontumour_sequencing_metadata[(nontumour_sequencing_metadata.fk_britroc_number == int(wildcards.patient_id))]

	# configure 'analysis_type' string
	if wildcards.analysis_type == 'panel_6_28':
		analysis_type = 'panel_6_28'
	elif wildcards.analysis_type == 'panel_28_only':
		analysis_type = 'antijoined_panel_28_6'
	elif wildcards.analysis_type == 'panel_28':
		analysis_type = 'panel_28'
	elif wildcards.analysis_type == 'TP53':
		analysis_type = 'union_of_tp53' 

	if bam_or_bai == 'bam':
		bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.nonoverlapping_id_place_holder.analysis_type_place_holder.bam', zip, SLX_ID=test_sample_metadata['fk_slx'], barcodes=test_sample_metadata['fk_barcode'], flowcell=test_sample_metadata['flowcell'], lane=test_sample_metadata['lane'])
	elif bam_or_bai == 'bai':
		bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.nonoverlapping_id_place_holder.analysis_type_place_holder.bai', zip, SLX_ID=test_sample_metadata['fk_slx'], barcodes=test_sample_metadata['fk_barcode'], flowcell=test_sample_metadata['flowcell'], lane=test_sample_metadata['lane'])

	bam_files = []

	for bam_file_name in bam_files_tmp:
		new_bam_name = bam_file_name.replace('nonoverlapping_id_place_holder', wildcards.nonoverlapping_id)
		new_bam_name = new_bam_name.replace('analysis_type_place_holder', analysis_type)
		bam_files.append(new_bam_name)

	return(bam_files)

def cleaned_normal_bams(wildcards):

	metadata = germline_metadata[(germline_metadata.fk_britroc_number == int(wildcards.patient_id))]

	# configure 'analysis_type' string
	if wildcards.analysis_type == 'panel_6_28':
		analysis_type = 'panel_6_28'
	elif wildcards.analysis_type == 'panel_28_only':
		analysis_type = 'antijoined_panel_28_6' 

	bam_files_tmp = expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.nonoverlapping_id_place_holder.analysis_type_place_holder.bam', zip, SLX_ID=metadata['fk_slx'], barcodes=metadata['fk_barcode'], flowcell=metadata['flowcell'], lane=metadata['lane'])

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

def get_nonoverlapping_id_list(wildcards):
	if wildcards.analysis_type == 'panel_6_28':
		nonoverlapping_id_list = [1,2,3,4]
	elif wildcards.analysis_type == 'panel_28_only':
		nonoverlapping_id_list = [1,2,3]
	return(nonoverlapping_id_list)

def get_gene_set_analysed(wildcards):
	if wildcards.analysis_type == 'panel_6_28':
		return(['TP53','BRCA1','BRCA2','FANCM','BARD1','RAD51B','RAD51C','RAD51D','BRIP1','PALB2'])
	elif wildcards.analysis_type == 'panel_28_only':
		return(['TP53','NRAS','PIK3CA','CTNNB1','EGFR','BRAF','PTEN','KRAS','RB1','CDK12','NF1'])
	elif wildcards.analysis_type == 'panel_28':
		return(["BRCA1","BRCA2","RAD51C","RAD51D","RAD51B","BRIP1","FANCM","PALB2","BARD1","CDK12","EGFR","PTEN","TP53","KRAS","BRAF","PIK3CA","CTNNB1","NF1","RB1","NRAS"])

def get_relevant_patient_list(wildcards):
	if wildcards.analysis_type=='panel_6_28':
		return(matched_somatic_patients)
	elif wildcards.analysis_type=='panel_28_only':
		return(matched_somatic_patients_panel_28_only)

def get_relevant_bed_file(wildcards):
	if wildcards.analysis_type == 'panel_6_28':
		output_file='resources/panel_6_28.nonoverlapping.targets.{}.bed'.format(wildcards.nonoverlapping_id)
	elif wildcards.analysis_type == 'panel_28_only':
		output_file='resources/nonoverlapping.antijoined_panel_28_6.targets.{}.bed'.format(wildcards.nonoverlapping_id)
	elif wildcards.analysis_type == 'panel_28':
		output_file='resources/nonoverlapping.panel_28.targets.{}.bed'.format(wildcards.nonoverlapping_id)	
	return(output_file)
