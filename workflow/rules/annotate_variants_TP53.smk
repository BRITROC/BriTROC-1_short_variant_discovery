rule vep_octopus:
	input:
		'results/tumour_sample_vcfs_octopus/tp53_final/{sample}.filtered5.vcf'
	output: 'results/tumour_sample_vcfs_octopus/tp53_final/{sample}.filtered.vep.vcf'
	conda: '../../config/vep.yaml'
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

#rule link_vep_filtering_script:
#	input: '../strelka/strelka_binary/collate_and_filter_vep_files.R'
#	output: 'workflow/scripts/collate_and_filter_vep_files.R'
#	shell: 'ln -s {input} {output}'

#rule collate_and_filter_archival_vep_files:
#	input: expand('results/tumour_sample_vcfs/{sample}.filtered.vep.vcf', sample=paired_archival_samples)
#	output: 'results/filtered_archival_vep_calls.tsv'
#	script: '../scripts/collate_and_filter_vep_files.R'

#rule collate_and_filter_archival_vardict_vep_files:
#	input: expand('results/tumour_sample_vcfs_vardict/{sample}.filtered.vep.vcf', sample=paired_archival_samples)
#	output: 'results/filtered_archival_vep_calls_vardict.tsv'
#	script: '../scripts/collate_and_filter_vep_files.R'

rule collate_and_filter_archival_octopus_vep_files:
	input: expand('results/tumour_sample_vcfs_octopus/{sample}.filtered.vep.vcf', sample=matched_archival_samples)
	output: 'results/final_tp53/filtered_archival_vep_calls_octopus_non_tp53.tsv'
	script: '../scripts/annotate_variants_TP53/collate_and_filter_vep_files_non_tp53.R'

# maybe we should perform joint-variant calling on all archival samples belonging to the same patient?

rule collate_and_filter_relapse_octopus_vep_files:
	input: expand('results/tumour_sample_vcfs_octopus/{sample}.filtered.vep.vcf', sample=somatic_tp53_samples)
	output: 'results/final_tp53/filtered_relapse_vep_calls_octopus.tsv'
	script: '../scripts/annotate_variants_TP53/collate_and_filter_vep_files.R'

rule collate_allele_fraction_data:
	input: expand('results/tumour_sample_vcfs_octopus/{sample}.filtered2.vcf', sample=somatic_tp53_samples)
	output: 'results/final_TP53/tp53_collated_MAFs.tsv'
	script: '../scripts/annotate_variants_TP53/extract_info_from_vcf.R'

rule add_MAFs_to_TP53_variant_table:
	input:
		filtered_archival_TP53_variants='results/final_tp53/filtered_archival_vep_calls_octopus.tsv',
		filtered_relapse_TP53_variants='results/final_tp53/filtered_relapse_vep_calls_octopus.tsv',
		TP53_variant_MAFs='results/final_tp53/tp53_collated_MAFs.tsv'
	output: 'results/final_tp53/filtered_TP53_variants_with_MAFs.tsv'
	script: '../scripts/annotate_variants_TP53/get_tp53_table.R'

rule classify_clonality_of_TP53_variants:
	input:
		filtered_TP53_variants_with_MAFs=rules.add_MAFs_to_TP53_variant_table.output.TP53_variant_MAFs,
	output: TP53_variants_classified_by_clonality='results/final_tp53/TP53_variants_with_clonality_classifications.tsv'
	script: '../scripts/annotate_variants_TP53/get_TP53_clonality_classifications.R'

#rule collate_and_filter_relapse_vep_files:
#	input: expand('results/tumour_sample_vcfs/{sample}.filtered.vep.vcf', sample=paired_relapse_samples)
#	output: 'results/filtered_relapse_vep_calls.tsv'
#	script: '../scripts/collate_and_filter_vep_files.R'

#rule get_shared_variants:
#	input: expand('results/tumour_sample_vcfs/{sample}.filtered.vcf', sample=paired_somatic_samples)
#	output: 'results/shared_calls.tsv'
#	script: '../scripts/get_shared_calls.R'
