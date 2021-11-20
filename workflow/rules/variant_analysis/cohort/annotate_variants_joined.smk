rule vep_octopus:
	input: rules.concat_vcfs.output
	output: 'results/variant_analysis/cohort/{patient_id}.filtered.vep.vcf'
	conda: '../../../../config/vep.yaml'
	shell: 'ensembl-vep/vep \
			-i {input} \
			-o {output} \
			--cache \
			--offline \
			--format vcf \
			--dir /Users/bradle02/.vep/ \
			--force_overwrite \
			--fasta /Users/bradle02/.vep/homo_sapiens/103_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
			--check_existing \
			--everything \
			--no_escape \
			-a GRCh37 \
			--port 3337'

rule ensure_tech_rep_genotypes_match:
	input: 
		all_patient_vcf=rules.concat_vcfs.output,
		metadata_sheet='config/all_tumour_metadata.tsv'
	output: 
		tumour_samples_union='results/variant_analysis/cohort/{patient_id}.filtered3.vcf',
		library_MAFs='results/variant_analysis/cohort/{patient_id}.library_MAFs.vcf',
		library_depths='results/variant_analysis/cohort/{patient_id}.library_depths.vcf',
		sample_genotypes='results/variant_analysis/cohort/{patient_id}.sample_genotypes.vcf'
	script: '../../../scripts/annotate_variants_joined/view_square_vcfs_cohort.R'

rule collate_and_filter_tumour_type_specific_vcf_files:
	input: lambda wildcards: expand('results/variant_analysis/cohort/{patient_id}.filtered3.vcf', patient_id=all_tumour_sample_patients)
	output: 'results/variant_analysis/cohort/collated/filtered3_joined.tsv'
	script: '../../../scripts/annotate_variants_joined/filtered4_files_joined.R'

rule collate_and_filter_octopus_vep_files:
	input: 
		vep_files= lambda wildcards: expand('results/variant_analysis/cohort/{patient_id}.filtered.vep.vcf',patient_id=all_tumour_sample_patients),
		vcf_file=rules.collate_and_filter_tumour_type_specific_vcf_files.output
	output: 'results/variant_analysis/cohort/collated/filtered_vep_calls_octopus_joined.tsv'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'

rule ensure_tech_rep_genotypes_match_with_tumour_type:
	input: 
		all_patient_vcf=rules.concat_vcfs.output,
		metadata_sheet='config/all_tumour_metadata.tsv'
	output: 
		tumour_samples_union='results/variant_analysis/cohort/both/{patient_id}.filtered3.vcf',
		archival_samples='results/variant_analysis/cohort/both/{patient_id}.archival.filtered3.vcf',
		relapse_samples='results/variant_analysis/cohort/both/{patient_id}.relapse.filtered3.vcf',
		library_MAFs='results/variant_analysis/cohort/both/{patient_id}.library_MAFs.vcf',
		library_depths='results/variant_analysis/cohort/both/{patient_id}.library_depths.vcf',
		sample_genotypes='results/variant_analysis/cohort/both/{patient_id}.sample_genotypes.vcf'
	script: '../../../scripts/annotate_variants_joined/view_square_vcfs_no_germline.R'

rule collate_and_filter_tumour_type_specific_vcf_files_with_tumour_type:
	input: lambda wildcards: expand('results/variant_analysis/cohort/both/{patient_id}.{tumour_type}.filtered3.vcf', patient_id=all_patients_with_tumour_samples_of_both_types, tumour_type=wildcards.tumour_type)
	output: 'results/variant_analysis/cohort/collated/{tumour_type}_filtered3_joined.tsv'
	script: '../../../scripts/annotate_variants_joined/filtered4_files_joined.R'

rule collate_and_filter_octopus_vep_files_with_tumour_type:
	input: 
		vep_files= lambda wildcards: expand('results/variant_analysis/cohort/{patient_id}.filtered.vep.vcf',patient_id=all_patients_with_tumour_samples_of_both_types),
		vcf_file=rules.collate_and_filter_tumour_type_specific_vcf_files_with_tumour_type.output
	output: 'results/variant_analysis/cohort/collated/filtered_{tumour_type}_vep_calls_octopus_joined.tsv'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'
