#rule vep_octopus_targeted:
#	input: rules.concat_vcfs_targeted.output
#	output: 'results/variant_analysis/unmatched/paired/{patient_id}.filtered.vep.targeted.vcf'
#	conda: '../../../../config/vep.yaml'
#	shell: 'ensembl-vep/vep \
#			-i {input} \
		#	-o {output} \
		#	--cache \
		#	--offline \
		#	--format vcf \
		#	--dir /Users/bradle02/.vep/ \
		#	--force_overwrite \
		#	--fasta /Users/bradle02/.vep/homo_sapiens/103_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
		#	--check_existing \
		#	--everything \
		#	--no_escape \
		#	-a GRCh37 \
		#	--port 3337'

rule ensure_tech_rep_genotypes_match_with_tumour_type_targeted:
	input: 
		combined_vcfs=rules.concat_vcfs_targeted.output,
		tumour_metadata='config/tumour_metadata_with_one_of_both_types.tsv'
	params:
		variant_quality_score_threshold=500,
		C_to_G_maf_threshold=0.23,
		C_to_G_maf_diff_threshold=0.30,
		includes_germline_sample_column=False,
		includes_tumour_type_analysis=True,
		includes_germline_variants=True
	output: 
		tumour_samples_union='results/variant_analysis/unmatched/paired/{patient_id}.filtered3.targeted.vcf',
		archival_samples='results/variant_analysis/unmatched/paired/{patient_id}.archival.filtered3.targeted.vcf',
		relapse_samples='results/variant_analysis/unmatched/paired/{patient_id}.relapse.filtered3.targeted.vcf',
		library_MAFs='results/variant_analysis/unmatched/paired/{patient_id}.library_MAFs.targeted.vcf',
		library_depths='results/variant_analysis/unmatched/paired/{patient_id}.library_depths.targeted.vcf',
		sample_genotypes='results/variant_analysis/unmatched/paired/{patient_id}.sample_genotypes.targeted.vcf'
	script: '../../../scripts/annotate_variants_joined/view_square_vcfs.R'

rule collate_and_filter_tumour_type_specific_vcf_files_with_tumour_type_targeted:
	input: lambda wildcards: expand('results/variant_analysis/unmatched/paired/{patient_id}.{tumour_type}.filtered3.targeted.vcf', patient_id=all_patients_with_tumour_samples_of_both_types, tumour_type=wildcards.tumour_type)
	output: 'results/variant_analysis/unmatched/paired/collated/{tumour_type}_filtered3_joined.targeted.tsv'
	script: '../../../scripts/annotate_variants_joined/filtered4_files_joined.R'

rule MTBP_filter_curated_results_tumour_type_specific_targeted:
	input: rules.collate_and_filter_tumour_type_specific_vcf_files_with_tumour_type_targeted.output
	output: 'results/variant_analysis/unmatched/paired/collated/{tumour_type}_filtered3_joined_MTBP_filtered.targeted.tsv'
	script: '../../../scripts/annotate_variants_joined/apply_MTBP_filter_unmatched_and_unpaired.R'

#rule collate_and_filter_octopus_vep_files_with_tumour_type_targeted:
#	input: 
#		vep_files= lambda wildcards: expand('results/variant_analysis/unmatched/paired/{patient_id}.filtered.vep.targeted.vcf',patient_id=all_patients_with_tumour_samples_of_both_types),
#		vcf_file=rules.collate_and_filter_tumour_type_specific_vcf_files_with_tumour_type_targeted.output
#	output: 'results/variant_analysis/unmatched/paired/collated/filtered_{tumour_type}_vep_calls_octopus_joined.targeted.tsv'
#	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'

rule collate_and_filter_octopus_vep_files_with_tumour_type_targeted:
	input:
		vep_files=lambda wildcards: expand('results/variant_analysis/unmatched/{patient_id}.filtered.vep.vcf',patient_id=all_tumour_sample_patients),
		vcf_file=rules.MTBP_filter_curated_results_tumour_type_specific_targeted.output
	output:
		vep_output='results/variant_analysis/unmatched/paired/collated/{tumour_type}_filtered_vep_calls_octopus_joined.targeted.tsv',
		vep_reduced='results/variant_analysis/unmatched/paired/collated/{tumour_type}_BriTROC-1_unmatched_and_paired_variants.targeted.tsv',
		vcf_output='results/variant_analysis/unmatched/paired/collated/{tumour_type}_filtered_calls_octopus_joined.targeted.vcf'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'
