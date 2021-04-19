def get_tumour_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = all_somatic_metadata[(all_somatic_metadata.fk_sample == wildcards.sample)].head(n=2)

	SLX = test_sample_metadata.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	barcodes = test_sample_metadata.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()
	lane = test_sample_metadata.set_index('lane', drop=False)
	lane = lane.index.unique().tolist()
	flowcell = test_sample_metadata.set_index('fk_run', drop=False)
	flowcell = flowcell.index.unique().tolist()
	flowcell = flowcell[0].split('_')

	bam_file_1 = '../SLX/{}/bam/cleaned_bams/{}.{}.{}.s_{}.amplicon.bam'.format(SLX[0],  ''.join(SLX), barcodes[0], flowcell[-1], lane[0])
	bam_file_2 = '../SLX/{}/bam/cleaned_bams/{}.{}.{}.s_{}.amplicon.bam'.format(SLX[0],  ''.join(SLX), barcodes[1], flowcell[-1], lane[0])

	return([bam_file_1, bam_file_2])

def get_normal_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = all_somatic_metadata[(all_somatic_metadata.fk_sample == wildcards.sample)].head(n=2)
	britroc_number = test_sample_metadata.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	patient_metadata = all_somatic_metadata[all_somatic_metadata['fk_britroc_number'] == britroc_number[0]]

	normal_metadata = germline_metadata[germline_metadata['fk_britroc_number'] == britroc_number[0]].head(n=2)

	SLX = normal_metadata.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()

	barcodes = normal_metadata.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()

	lane = normal_metadata.set_index('lane', drop=False)
	lane = lane.index.unique().tolist()
	
	flowcell = normal_metadata.set_index('fk_run', drop=False)
	flowcell = flowcell.index.unique().tolist()
	flowcell = flowcell[0].split('_')

	if SLX[0] == 'SLX-9629':

		bam_file_1 = '../SLX/{}/samplebam/cleaned_samplebams/{}.{}.amplicon.bam'.format(SLX[0], SLX[0], barcodes[0])
		return([bam_file_1])

	else:
		bam_file_1 = '../SLX/{}/bam/cleaned_bams/{}.{}.{}.s_{}.amplicon.bam'.format(SLX[0],  ''.join(SLX), barcodes[0], flowcell[-1], lane[0])
		bam_file_2 = '../SLX/{}/bam/cleaned_bams/{}.{}.{}.s_{}.amplicon.bam'.format(SLX[0],  ''.join(SLX), barcodes[1], flowcell[-1], lane[0])

		return([bam_file_1, bam_file_2])

def get_normal_sample_names(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = all_somatic_metadata[(all_somatic_metadata.fk_sample == wildcards.sample)].head(n=2)
	britroc_number = test_sample_metadata.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	patient_metadata = all_somatic_metadata[all_somatic_metadata['fk_britroc_number'] == britroc_number[0]]
	
	normal_metadata = germline_metadata[germline_metadata['fk_britroc_number'] == britroc_number[0]].head(n=2)

	SLX = normal_metadata.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	
	samples = normal_metadata.set_index('fk_sample', drop=False)
	samples = samples.index.unique().tolist()	

	barcodes = normal_metadata.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()

	if SLX[0] == 'SLX-9629':
		samples[0] = samples[0].replace('-','')
	elif SLX[0] in ['SLX-13630','SLX-11347','SLX-11111','SLX-11109','SLX-9856']:
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

rule convert_interval_list_to_bed_format:
	input: 
		targets='resources/{amplicon_panel_id}_amplicons.targets.interval_list',
		amplicons='resources/{amplicon_panel_id}_amplicons.amplicons.interval_list'
	output: 'resources/{amplicon_panel_id}.bed'
	script: '../scripts/convert_picard_interval_to_bed.R'

rule convert_bed_to_oct_format:
	input: rules.convert_interval_list_to_bed_format.output
	output: 'resources/{amplicon_panel_id}.oct'
	script: '../scripts/convert_bed_to_octopus.R'

rule octopus:
	input:
		reference_genome=config['reference_genome'],
		interval_file='resources/intersected_panel_6_28.oct',
		tumour_bams=get_tumour_bam_files,
		normal_bams=get_normal_bam_files
	output: 
		tumour_vcf='results/tumour_sample_vcfs_octopus/{sample}.vcf',
	threads: 4
	params:
		normal_sample_identifier=get_normal_sample_names
	run:
			shell ('../octopus/bin/octopus \
				-C cancer \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--threads \
				-w temp/ \
				-I {input.tumour_bams[0]} {input.tumour_bams[1]} {input.normal_bams[0]} \
				-N {params.normal_sample_identifier[0]} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}')

rule filter_octopus_calls:
	input:
		reference_genome=config['reference_genome'],
		interval_file='resources/intersected_panel_6_28.oct',
		tumour_bams=get_tumour_bam_files,
		normal_bams=get_normal_bam_files,
		vcf_file='results/tumour_sample_vcfs_octopus/{sample}.filtered.vcf'
	output: 
		tumour_vcf='results/tumour_sample_vcfs_octopus/{sample}.filtered2.vcf',
	threads: 4
	params:
		normal_sample_identifier=get_normal_sample_names
	run:
			shell ('../octopus/bin/octopus \
				-C cancer \
				--allow-octopus-duplicates \
				--threads \
				-w temp/ \
				--filter-vcf {input.vcf_file} \
				--somatics-only \
				--somatic-filter-expression "QUAL < 20 | AF < 0.05 | DP < 10" \
				-I {input.tumour_bams[0]} {input.tumour_bams[1]} {input.normal_bams[0]} \
				-N {params.normal_sample_identifier[0]} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}')

rule filter_octopus_calls2:
	input:
		reference_genome=config['reference_genome'],
		interval_file='resources/intersected_panel_6_28.oct',
		tumour_bams=get_tumour_bam_files,
		normal_bams=get_normal_bam_files,
		vcf_file='results/tumour_sample_vcfs_octopus/{sample}.filtered3.vcf'
	output: 
		tumour_vcf='results/tumour_sample_vcfs_octopus/{sample}.filtered4.vcf',
	threads: 4
	params:
		normal_sample_identifier=get_normal_sample_names
	run:
			shell ('../octopus/bin/octopus \
				-C cancer \
				--allow-octopus-duplicates \
				--threads \
				-w temp/ \
				--filter-vcf {input.vcf_file} \
				--somatic-forest resources/somatic.v0.7.2.forest \
				--somatics-only \
				-I {input.tumour_bams[0]} {input.tumour_bams[1]} {input.normal_bams[0]} \
				-N {params.normal_sample_identifier[0]} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}')


# this is a fudge so that the 'match genotype' script works
rule remove_normal_sample:
	input: 'results/tumour_sample_vcfs_octopus/{sample}.filtered5.vcf'
	output: 'results/tumour_sample_vcfs_octopus/{sample}.filtered6.vcf'
	params: normal_sample_identifier=get_normal_sample_names
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -s "^{params.normal_sample_identifier[0]}" {input} > {output}'

#rule symlink_results:
#	input: rules.octopus.output.tumour_vcf
#	output: 'results/tumour_sample_vcfs_octopus/{sample}.vcf'
#	shell: 'ln {input} {output}'

rule analyse_tp53:
	input: rules.octopus.output.tumour_vcf
	output: 'results/tumour_sample_vcfs_octopus/{sample}_tp53.txt'
	script: '../scripts/process_oct_output.R'
