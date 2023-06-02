rule vep_octopus_germline:
	input: rules.ensure_genotypes_match.output
	output: 'results/variant_analysis/germline/octopus_unmatched/{analysis_type}/merged/{patient_id}.vep.vcf'
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

rule collate_and_filter_octopus_germline_vcfs:
	input: 
		vcf_files= lambda wildcards: expand('results/variant_analysis/germline/octopus_unmatched/{analysis_type}/merged/{patient_id}.vcf',patient_id=patients_with_nontumour_sample_sequencing, analysis_type=wildcards.analysis_type)
	output: 
		vcf_full='results/variant_analysis/germline/octopus_unmatched/{analysis_type}/merged/collated/octopus_unmatched_germline_collated.tsv',
		vcf_MTBP='results/variant_analysis/germline/octopus_unmatched/{analysis_type}/merged/collated/octopus_unmatched_germline_collated_MTBP.vcf'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vcf_files_joined_octopus_germline.R'

rule MTBP_curate_germline_octopus_vcfs:
	input: rules.collate_and_filter_octopus_germline_vcfs.output.vcf_full
	output: 'results/variant_analysis/germline/octopus_unmatched/{analysis_type}/merged/collated/octopus_unmatched_germline_collated_final.vcf'
	script: '../../../scripts/MTBP_curate_germline_variants_octopus.R'

rule collate_and_filter_octopus_vep_files_germline:
	input: 
		vep_files= lambda wildcards: expand('results/variant_analysis/germline/octopus_unmatched/{analysis_type}/merged/{patient_id}.vep.vcf',patient_id=patients_with_nontumour_sample_sequencing, analysis_type=wildcards.analysis_type),
		vcf_file=rules.MTBP_curate_germline_octopus_vcfs.output
	output: 
		vep_output='results/variant_analysis/germline/octopus_unmatched/{analysis_type}/merged/collated/filtered_vep_calls_octopus_joined.tsv',
		vep_reduced='results/variant_analysis/germline/octopus_unmatched/{analysis_type}/merged/collated/BriTROC-1_matched_and_unpaired_somatic_variants.tsv',
		vcf_output='results/variant_analysis/germline/octopus_matched/{analysis_type}/merged/collated/filtered_calls_octopus_joined.vcf'
	script: '../../../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'
