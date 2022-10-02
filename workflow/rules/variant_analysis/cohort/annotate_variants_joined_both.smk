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

rule ensure_tech_rep_genotypes_match_with_tumour_type:
	input: 
		combined_vcfs='results/variant_analysis/cohort/{patient_id}.filtered2.vcf',
		tumour_metadata='config/tumour_metadata_with_one_of_both_types.tsv'
	params:
		variant_quality_score_threshold=500,
		C_to_G_maf_threshold=0.23,
		C_to_G_maf_diff_threshold=0.30,
		includes_germline_sample_column=False,
		includes_tumour_type_analysis=True,
		includes_germline_variants=True
	output: 
		tumour_samples_union='results/variant_analysis/cohort/both/{patient_id}.filtered3.vcf',
		archival_samples='results/variant_analysis/cohort/both/{patient_id}.archival.filtered3.vcf',
		relapse_samples='results/variant_analysis/cohort/both/{patient_id}.relapse.filtered3.vcf',
		library_MAFs='results/variant_analysis/cohort/both/{patient_id}.library_MAFs.vcf',
		library_depths='results/variant_analysis/cohort/both/{patient_id}.library_depths.vcf',
		sample_genotypes='results/variant_analysis/cohort/both/{patient_id}.sample_genotypes.vcf'
	script: '../../../scripts/annotate_variants_joined/view_square_vcfs.R'

rule get_interval_file_for_targeted_calling:
	input: rules.ensure_tech_rep_genotypes_match_with_tumour_type.output.tumour_samples_union
	output: 'results/variant_analysis/cohort/both/{patient_id}.interval_file.txt'
	script: '../../../scripts/annotate_variants_joined/get_intervals_for_targeted_calling.R'

rule octopus_targeted_calling:
	input:
		reference_genome=config['reference_genome'],
		interval_file=rules.get_interval_file_for_targeted_calling.output, 
		tumour_bams=get_tumour_bam_files,
		tumour_bam_indexes=get_tumour_bam_index_files
	output: 
		tumour_vcf='results/variant_analysis/cohort/{patient_id}.{nonoverlapping_id}.targeted.vcf',
	threads: 16
	wildcard_constraints:
		nonoverlapping_id='[1-9]'
	shell: 'octopus \
				-C cancer \
				--disable-downsampling \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--forest resources/germline.v0.7.2.forest \
				--somatic-forest resources/somatic.v0.7.2.forest \
				--max-somatic-haplotypes 2 \
				--annotations SB SD AF AD FRF \
				--filter-expression "QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.50 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--somatic-filter-expression "QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.90 | SD > 0.90 | FRF > 0.5 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | AD < 1 | AF < 0.0001" \
				--min-expected-somatic-frequency 0.03 \
				--min-credible-somatic-frequency 0.01 \
				--threads \
				-w temp/ \
				-I {input.tumour_bams} \
				--regions-file {input.interval_file} \
				--output {output.tumour_vcf} \
				-R {input.reference_genome}'

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
