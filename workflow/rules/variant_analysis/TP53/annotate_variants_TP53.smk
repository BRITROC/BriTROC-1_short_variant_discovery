rule vep_octopus_tp53:
	input: rules.concat_vcfs_sample_level.output
	output: 'results/variant_analysis/{analysis_type}/{sample}.filtered.vep.vcf'
	conda: '../../../../config/vep.yaml'
	shell: 'ensembl-vep/vep \
			-i {input} \
			-o {output} \
			--cache \
			--offline \
			--format vcf \
			--dir vep_cache/ \
			--force_overwrite \
			--fasta /Users/bradle02/.vep/homo_sapiens/103_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
			--check_existing \
			--everything \
			--no_escape \
			-a GRCh37 \
			--port 3337'

rule collate_and_filter_octopus_vep_files_tp53:
	input: 
		filtered_vep_files=expand('results/variant_analysis/TP53/{sample}.filtered.vep.vcf', sample=TP53_sequenced_DNA_samples),
		britroc1_DNA_samples='britroc1_db/database_text_file_output/sample.tsv'
	output: 'results/variant_analysis/TP53/collated/filtered_{tumour_type}_vep_calls_octopus.tsv'
	script: '../../../scripts/annotate_variants_TP53/collate_and_filter_vep_files.R'

rule collate_allele_fraction_data:
	input: expand('results/variant_analysis/TP53/{sample}.vcf', sample=TP53_sequenced_DNA_samples)
	output: 'results/variant_analysis/TP53/collated/tp53_collated_MAFs.tsv'
	script: '../../../scripts/annotate_variants_TP53/extract_info_from_vcf.R'

rule add_MAFs_to_TP53_variant_table:
	input:
		filtered_archival_TP53_variants='results/variant_analysis/TP53/collated/filtered_archival_vep_calls_octopus.tsv',
		filtered_relapse_TP53_variants='results/variant_analysis/TP53/collated/filtered_relapse_vep_calls_octopus.tsv',
		TP53_variant_MAFs='results/variant_analysis/TP53/collated/tp53_collated_MAFs.tsv'
	output: 'results/variant_analysis/TP53/collated/filtered_TP53_variants_with_MAFs.tsv'
	script: '../../../scripts/annotate_variants_TP53/get_tp53_table.R'

rule classify_clonality_of_TP53_variants:
	input:
		filtered_TP53_variants_with_MAFs=rules.add_MAFs_to_TP53_variant_table.output,
	output: TP53_variants_classified_by_clonality='results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv'
	script: '../../../scripts/annotate_variants_TP53/get_TP53_clonality_classifications.R'
