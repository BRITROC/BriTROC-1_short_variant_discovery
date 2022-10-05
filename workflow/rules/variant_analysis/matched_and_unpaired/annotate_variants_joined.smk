rule vep_octopus:
	input: rules.concat_vcfs.output
	output: 'results/variant_analysis/matched_and_unpaired/{analysis_type}/{patient_id}.filtered.vep.vcf'
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
		combined_vcfs=rules.concat_vcfs.output,
		germline_metadata='config/germline_metadata_panel_matched_and_unpaired.tsv',
		tumour_metadata='config/somatic_metadata_panel_matched_and_unpaired.tsv'
	params:
		variant_quality_score_threshold=500,
		C_to_G_maf_threshold=0.23,
		C_to_G_maf_diff_threshold=0.30,
		includes_germline_sample_column=True,
		includes_tumour_type_analysis=False,
		includes_germline_variants=False
	output: 
		tumour_samples_union='results/variant_analysis/matched_and_unpaired/{analysis_type}/{patient_id}.filtered3.vcf',
		library_MAFs='results/variant_analysis/matched_and_unpaired/{analysis_type}/{patient_id}.library_MAFs.vcf',
		library_depths='results/variant_analysis/matched_and_unpaired/{analysis_type}/{patient_id}.library_depths.vcf',
		sample_genotypes='results/variant_analysis/matched_and_unpaired/{analysis_type}/{patient_id}.sample_genotypes.vcf'
	script: '../../../scripts/annotate_variants_joined/view_square_vcfs.R'

rule collate_and_filter_tumour_type_specific_vcf_files:
	input: lambda wildcards: expand('results/variant_analysis/matched_and_unpaired/{analysis_type}/{patient_id}.filtered3.vcf', patient_id=matched_and_unpaired_somatic_metadata_patients, analysis_type=wildcards.analysis_type)
	output: 'results/variant_analysis/matched_and_unpaired/{analysis_type}/collated/filtered3_joined.tsv'
	script: '../../../scripts/annotate_variants_joined/filtered4_files_joined.R'

rule collate_and_filter_octopus_vep_files:
	input: 
		vep_files= lambda wildcards: expand('results/variant_analysis/matched_and_unpaired/{analysis_type}/{patient_id}.filtered.vep.vcf',patient_id=matched_and_unpaired_somatic_metadata_patients, analysis_type=wildcards.analysis_type),
		vcf_file=rules.collate_and_filter_tumour_type_specific_vcf_files.output
	output: 'results/variant_analysis/matched_and_unpaired/{analysis_type}/collated/filtered_vep_calls_octopus_joined.tsv'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'
