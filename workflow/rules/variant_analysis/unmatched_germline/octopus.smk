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

#def return_group_of_vcf_files(wildcards):
#	vcf_metadata = pandas.read_table("results/variant_analysis/germline/vcf_table.tsv").set_index('vcf_name', drop=False)
#	vcf_metadata = vcf_metadata[(vcf_metadata.group == int(wildcards.group_number))]
#
#	vcf_files = vcf_metadata.index.unique().tolist()
#
#	return(vcf_files)

#def return_group_of_bam_files(wildcards):
#	vcf_metadata = pandas.read_table("results/variant_analysis/germline/vcf_table.tsv").set_index('bam_name', drop=False)
#	vcf_metadata = vcf_metadata[(vcf_metadata.group == int(wildcards.group_number))]
#	bam_files = vcf_metadata.index.unique().tolist()
#
#	return(bam_files)

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

rule octopus_germline:
	input:
		interval_file=rules.convert_bed6_to_oct_format.output,
		normal_bams=cleaned_normal_bams,
		reference_genome=config['reference_genome']
	threads: 4
	output: tumour_vcf=protected('results/variant_analysis/germline/octopus_unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.vcf')
	shell: '../octopus/bin/octopus \
			-R {input.reference_genome} \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--disable-downsampling \
				--annotations SB SD AF AD FRF \
				--filter-expression "QUAL < 100 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.25 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--threads \
				--regions-file {input.interval_file} \
				-I {input.normal_bams} -o {output.tumour_vcf}'

rule filter_octopus_raw_calls_germline:
	input: 
		filtered_vcf=rules.octopus_germline.output
	output: 'results/variant_analysis/germline/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.1.panel_6_28.passed.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

#rule octopus_germline_joint_calling:
#	input: single_library_vcfs=return_group_of_vcf_files,
#		reference_genome=config['reference_genome'],
#		normal_bams=return_group_of_bam_files
#	threads: 15
#	output: 'results/variant_analysis/germline/jointly_called.group_{group_number}.vcf'
#	shell: '../octopus/bin/octopus \
		#			-R {input.reference_genome} \
		#	--annotations AF AD \
		#	-I {input.normal_bams} \
		#	--disable-denovo-variant-discovery \
		#	--threads \
		#	-c {input.single_library_vcfs} \
		#	-o {output}'

#rule octopus_germline_joint_calling_2:
#	input: single_library_vcfs=expand('results/variant_analysis/germline/jointly_called.group_{group_number}.vcf', group_number=list(range(1,21,1))),
#		reference_genome=config['reference_genome'],
#		normal_bams=expand('../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.1.panel_6_28.bam', zip, SLX_ID=matched_germline_metadata['fk_slx'], barcodes=matched_germline_metadata['fk_barcode'], flowcell=matched_germline_metadata['flowcell'], lane=matched_germline_metadata['lane'])
#	threads: 15
#	output: 'results/variant_analysis/germline/jointly_called.vcf'
#	shell: '../octopus/bin/octopus \
		#			-R {input.reference_genome} \
		#	--annotations AF AD \
		#	-I {input.normal_bams} \
		#	--disable-denovo-variant-discovery \
		#	--threads \
		#	--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.25 | SB > 0.98 | BQ < 15 | DP < 21 | ADP < 1" \
		#	--max-genotype-combinations 100000 \
		#	-c {input.single_library_vcfs} \
		#	-o {output}'
