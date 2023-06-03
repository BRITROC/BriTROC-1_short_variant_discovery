# run deep variant on whole blood samples

def get_normal_bam_files(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = germline_metadata[(germline_metadata.fk_britroc_number == int(wildcards.patient_id))]
	britroc_number = test_sample_metadata.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

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

def get_normal_bam_files_container(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	test_sample_metadata = germline_metadata[(germline_metadata.fk_britroc_number == int(wildcards.patient_id))]
	britroc_number = test_sample_metadata.set_index('fk_britroc_number', drop=False)
	britroc_number = britroc_number.index.unique().tolist()

	# configure 'analysis_type' string
	if wildcards.analysis_type == 'panel_6_28':
		analysis_type = 'panel_6_28'
	elif wildcards.analysis_type == 'panel_28_only':
		analysis_type = 'antijoined_panel_28_6' 

	bam_files_tmp = expand('/bam_dir/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.nonoverlapping_id_place_holder.analysis_type_place_holder.bam', zip, SLX_ID=test_sample_metadata['fk_slx'], barcodes=test_sample_metadata['fk_barcode'], flowcell=test_sample_metadata['flowcell'], lane=test_sample_metadata['lane']) 
	bam_files = []

	for bam_file_name in bam_files_tmp:
		new_bam_name = bam_file_name.replace('nonoverlapping_id_place_holder', wildcards.nonoverlapping_id)
		new_bam_name = new_bam_name.replace('analysis_type_place_holder', analysis_type)
		bam_files.append(new_bam_name)

	return(bam_files)

def get_relevant_bed_file(wildcards):
	if wildcards.analysis_type == 'panel_6_28':
		output_file='resources/panel_6_28.nonoverlapping.targets.{}.bed'.format(wildcards.nonoverlapping_id)
	elif wildcards.analysis_type == 'panel_28_only':
		output_file='resources/nonoverlapping.antijoined_panel_28_6.amplicons.{}.bed'.format(wildcards.nonoverlapping_id)
	
	return(output_file)

#rule convert_bed6_to_oct_format:
#	input:  get_relevant_bed_file                                           
#	output: 'resources/{analysis_type}.{nonoverlapping_id}.targets.oct'
#	script: '../scripts/octopus_formatting/convert_bed6_to_octopus.R'

rule deepvariant:
	input:
		reference_genome=config['reference_genome'],
		interval_file='resources/panel_6_28.nonoverlapping.amplicons.{nonoverlapping_id}.bed', 
		normal_bams=get_normal_bam_files
	output: 
		vcf='results/variant_analysis/unmatched_germline/{analysis_type}/deepvariant/{patient_id}.{nonoverlapping_id}.vcf.gz',
		gvcf='results/variant_analysis/unmatched_germline/{analysis_type}/deepvariant/{patient_id}.{nonoverlapping_id}.g.vcf.gz'
	threads: 4
	params:
		container_input_bam_paths=get_normal_bam_files_container
	wildcard_constraints:
		nonoverlapping_id='[1-9]'
	container: 'workflow/rules/variant_analysis/unmatched_germline/deepvariant_1.4.0.sif'
	shell: '/opt/deepvariant/bin/run_deepvariant \
				--model_type=WES \
				--ref=/fasta_dir/hsa.GRCh37_g1kp2.fa \
				--reads={params.container_input_bam_paths} \
				--regions={input.interval_file} \
				--output_vcf={output.vcf} \
				--output_gvcf={output.gvcf} \
				--num_shards=4'


rule merge_deep_variant_results:
	input:
		reference_genome=config['reference_genome'],
		interval_file='resources/panel_6_28.nonoverlapping.amplicons.1.bed', 
		gvcf= lambda wildcards: expand('results/variant_analysis/unmatched_germline/{analysis_type}/deepvariant/{patient_id}.{nonoverlapping_id}.g.vcf.gz', analysis_type=wildcards.analysis_type, patient_id=germline_patients, nonoverlapping_id=wildcards.nonoverlapping_id)
	output: 
		gvcf='results/variant_analysis/unmatched_germline/{analysis_type}/deepvariant/all_patients.{nonoverlapping_id}.g.vcf.gz'
	threads: 4
	wildcard_constraints:
		nonoverlapping_id='[1-9]'
	container: 'workflow/rules/variant_analysis/unmatched_germline/glnexus_v1.3.1.sif'
	shell: '/usr/local/bin/glnexus_cli \
			--config DeepVariantWES \
			--bed {input.interval_file} \
			{input.gvcf} \
			| bcftools view - \
			| bgzip -c > {output.gvcf}'

rule index_deep_variant_compressed_vcf:
	input: rules.merge_deep_variant_results.output
	output: 'results/variant_analysis/unmatched_germline/{analysis_type}/deepvariant/all_patients.{nonoverlapping_id}.g.vcf.gz.csi'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools index {input}'

rule concat_deepvariant_vcfs:
	input:
		compressed_vcfs=lambda wildcards: expand('results/variant_analysis/unmatched_germline/{analysis_type}/deepvariant/all_patients.{nonoverlapping_id}.g.vcf.gz', analysis_type=wildcards.analysis_type, nonoverlapping_id=get_nonoverlapping_id_list(wildcards)),
		compressed_vcf_indexes=lambda wildcards: expand('results/variant_analysis/unmatched_germline/{analysis_type}/deepvariant/all_patients.{nonoverlapping_id}.g.vcf.gz.csi', analysis_type=wildcards.analysis_type, nonoverlapping_id=get_nonoverlapping_id_list(wildcards))
	output: protected('results/variant_analysis/unmatched_germline/{analysis_type}/deepvariant/all_patients.g.vcf.gz')
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools concat --allow-overlaps {input.compressed_vcfs} -O v -o {output}'

rule vep_deep_variant:
	input: rules.concat_deepvariant_vcfs.output
	output: 'results/variant_analysis/unmatched_germline/{analysis_type}/deepvariant/all_patients.vep.vcf'
	conda: '../../../../config/vep.yaml'
	shell: 'ensembl-vep/vep \
			-i {input} \
			-o {output} \
			--cache \
			--offline \
			--format vcf \
			--dir /Users/bradle02/.vep/ \
			--force_overwrite \
			--hgvsg \
			--fasta /Users/bradle02/.vep/homo_sapiens/103_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
			--check_existing \
			--everything \
			--no_escape \
			-a GRCh37 \
			--port 3337'

rule collate_and_filter_deep_variant_vep_files:
	input: 
		vep_files=rules.vep_deep_variant.output
	output: 'results/variant_analysis/unmatched_germline/{analysis_type}/deepvariant/all_patients.vep.filtered.vcf'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_deep_variant.R'

	#shell: 'singularity run \
	#			--bind /scratchb/jblab/bradle02/britroc1/tamseq/SLX/:/bam_dir \
	#			--bind /mnt/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh37_g1kp2/fasta/:/fasta_dir \
	#			workflow/rules/variant_analysis/unmatched_germline/deepvariant_1.4.0.sif \
	#			/opt/deepvariant/bin/run_deepvariant \
	#			--model_type=WES \
	#			--ref={input.reference_genome} \
	#			--reads={params.container_input_bam_paths} \
	#			--regions={input.interval_file} \
	#			--output_vcf={output.vcf} \
	#			--output_gvcf={output.gvcf} \
	#			--num_shards={threads}' 
