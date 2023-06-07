rule generate_vep_annotations_for_unmatched_variants:
	input: rules.concat_unmatched_vcfs_across_amplicon_groups.output
	output: 'results/variant_analysis/unmatched/{analysis_type}/{patient_id}.filtered.vep.vcf'
	conda: '../../../../config/vep.yaml'
	shell: 'ensembl-vep/vep \
			-i {input} \
			-o {output} \
			--cache \
			--offline \
			--format vcf \
			--dir vep_cache \
			--force_overwrite \
			--fasta /Users/bradle02/.vep/homo_sapiens/103_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
			--check_existing \
			--everything \
			--no_escape \
			-a GRCh37 \
			--port 3337'

rule ensure_tech_rep_genotypes_match_for_unmatched_variants:
	input:
		combined_vcfs=rules.concat_unmatched_vcfs_across_amplicon_groups.output,
		tumour_metadata='config/panel_28_tumour_amplicon_sequencing_metadata.tsv'
	params:
		variant_quality_score_threshold=500,
		C_to_G_maf_threshold=0.23,
		C_to_G_maf_diff_threshold=0.30,
		includes_germline_sample_column=False,
		includes_tumour_type_analysis=False,
		includes_germline_variants=True
	output: 
		tumour_samples_union='results/variant_analysis/unmatched/{analysis_type}/{patient_id}.filtered3.vcf',
		library_MAFs='results/variant_analysis/unmatched/{analysis_type}/{patient_id}.library_MAFs.vcf',
		library_depths='results/variant_analysis/unmatched/{analysis_type}/{patient_id}.library_depths.vcf',
		sample_genotypes='results/variant_analysis/unmatched/{analysis_type}/{patient_id}.sample_genotypes.vcf'
	script: '../../../scripts/annotate_variants_joined/view_square_vcfs.R'

rule collate_unmatched_variants_across_patients:
	input: lambda wildcards: expand('results/variant_analysis/unmatched/{analysis_type}/{patient_id}.filtered3.vcf', patient_id=patients_with_panel_28_tumour_sequencing, analysis_type=wildcards.analysis_type)
	output: 'results/variant_analysis/unmatched/{analysis_type}/collated/filtered3_joined.tsv'
	script: '../../../scripts/annotate_variants_joined/filtered4_files_joined.R'

rule reformat_vcf_for_MTBP:
	input: rules.collate_unmatched_variants_across_patients.output
	output: 'results/variant_analysis/unmatched/{analysis_type}/collated/filtered3_joined_MTBP_format.tsv'
	script: '../../../scripts/annotate_variants_joined/get_MTBP_format.R'

# TODO: Abstract MTBP annotations as a seperate data sheet
rule MTBP_filter_unmatched_variants:
	input: rules.collate_unmatched_variants_across_patients.output
	output: 'results/variant_analysis/unmatched/{analysis_type}/collated/filtered3_joined_MTBP_filtered.tsv'
	script: '../../../scripts/annotate_variants_joined/apply_MTBP_filter_unmatched_and_unpaired.R'

rule filter_unmatched_vcfs_using_vep_annotations:
	input: 
		vep_files= lambda wildcards: expand('results/variant_analysis/unmatched/{analysis_type}/{patient_id}.filtered.vep.vcf',patient_id=patients_with_panel_28_tumour_sequencing, analysis_type=wildcards.analysis_type),
		vcf_file=rules.MTBP_filter_unmatched_variants.output
	output: 
		vep_output='results/variant_analysis/unmatched/{analysis_type}/collated/filtered_vep_calls_octopus_joined.tsv',
		vep_reduced='results/variant_analysis/unmatched/{analysis_type}/collated/BriTROC-1_unmatched_and_unpaired_variants.tsv',
		vcf_output='results/variant_analysis/unmatched/{analysis_type}/collated/filtered_calls_octopus_joined.vcf'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'
