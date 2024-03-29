# copyright Thomas Bradley 2023 ('thomas.bradley@cruk.cam.ac.uk')

rule vep_octopus:
	input: 
		 variants=rules.concat_vcfs.output,
		 reference_genome=rules.decompress_reference_genome_for_use_with_vep.output
	output: 'results/variant_analysis/matched/{analysis_type}/{patient_id}.filtered.vep.vcf'
	conda: '../../../../config/vep.yaml'
	shell: 'ensembl-vep/vep \
			-i {input} \
			-o {output} \
			--cache \
			--offline \
			--format vcf \
			--dir vep_cache \
			--hgvsg \
			--force_overwrite \
			--fasta {input.reference_genome} \
			--check_existing \
			--everything \
			--no_escape \
			-a GRCh37 \
			--port 3337'

rule get_somatic_variants_only:
	input: rules.concat_vcfs.output
	output: 'results/variant_analysis/matched/{analysis_type}/{patient_id}.filtered2.somatics_only.vcf'
	shell: 'grep -e "SOMATIC" -e "#" {input} > {output}'

rule ensure_tech_rep_genotypes_match:
	input: 
		combined_vcfs=rules.get_somatic_variants_only.output,
		germline_metadata='config/nontumour_amplicon_sequencing_metadata.tsv',
		tumour_metadata='config/panel_28_tumour_amplicon_sequencing_metadata.tsv'
	params:
		variant_quality_score_threshold=500,
		C_to_G_maf_threshold=0.23,
		C_to_G_maf_diff_threshold=0.30,
		includes_germline_sample_column=True,
		includes_tumour_type_analysis=False,
		includes_germline_variants=False
	output: 
		tumour_samples_union='results/variant_analysis/matched/{analysis_type}/{patient_id}.filtered3.vcf',
		library_MAFs='results/variant_analysis/matched/{analysis_type}/{patient_id}.library_MAFs.vcf',
		library_depths='results/variant_analysis/matched/{analysis_type}/{patient_id}.library_depths.vcf',
		sample_genotypes='results/variant_analysis/matched/{analysis_type}/{patient_id}.sample_genotypes.vcf'
	script: '../../../scripts/annotate_variants_joined/view_square_vcfs.R'

rule collate_and_filter_vcf_files:
	input: lambda wildcards: expand('results/variant_analysis/matched/{analysis_type}/{patient_id}.filtered3.vcf', patient_id=patients_with_matched_and_unpaired_sequencing, analysis_type=wildcards.analysis_type)
	output: 'results/variant_analysis/matched/{analysis_type}/collated/filtered3_joined.tsv'
	script: '../../../scripts/annotate_variants_joined/filtered4_files_joined.R'

rule reformat_somatic_vcf_for_MTBP:
	input: rules.collate_and_filter_vcf_files.output
	output: 'results/variant_analysis/matched/{analysis_type}/collated/filtered3_joined_MTBP_format.tsv'
	script: '../../../scripts/annotate_variants_joined/get_MTBP_format.R'

rule MTBP_filter_curated_results:
	input: rules.collate_and_filter_vcf_files.output
	output: 'results/variant_analysis/matched/{analysis_type}/collated/filtered3_joined_MTBP_filtered.tsv'
	script: '../../../scripts/annotate_variants_joined/apply_MTBP_filter.R'

rule collate_and_filter_octopus_vep_files:
	input: 
		vep_files= lambda wildcards: expand('results/variant_analysis/matched/{analysis_type}/{patient_id}.filtered.vep.vcf',patient_id=patients_with_matched_and_unpaired_sequencing, analysis_type=wildcards.analysis_type),
		vcf_file=rules.MTBP_filter_curated_results.output
	output: 
		vep_output='results/variant_analysis/matched/{analysis_type}/collated/filtered_vep_calls_octopus_joined.tsv',
		vep_reduced='results/variant_analysis/matched/{analysis_type}/collated/BriTROC-1_matched_and_unpaired_somatic_variants.tsv',
		vcf_output='results/variant_analysis/matched/{analysis_type}/collated/filtered_calls_octopus_joined.vcf'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'
