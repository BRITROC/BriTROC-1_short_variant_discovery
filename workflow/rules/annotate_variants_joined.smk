rule vep_octopus:
	input:
		'results/tumour_sample_vcfs_octopus/{patient_id}.filtered2.vcf'
	output: 'results/tumour_sample_vcfs_octopus/{patient_id}.filtered.vep.vcf'
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

rule ensure_tech_rep_genotypes_match:
	input: 'results/tumour_sample_vcfs_octopus/{patient_id}.filtered2.vcf'
	output: 
		tumour_samples_union='results/tumour_sample_vcfs_octopus/{patient_id}.filtered3.vcf',
		archival_samples='results/tumour_sample_vcfs_octopus/{patient_id}.archival.filtered3.vcf',
		relapse_samples='results/tumour_sample_vcfs_octopus/{patient_id}.relapse.filtered3.vcf',
		library_MAFs='results/tumour_sample_vcfs_octopus/{patient_id}.library_MAFs.vcf',
		library_depths='results/tumour_sample_vcfs_octopus/{patient_id}.library_depths.vcf',
		sample_genotypes='results/tumour_sample_vcfs_octopus/{patient_id}.sample_genotypes.vcf'
	script: '../scripts/view_square_vcfs.R'

rule collate_and_filter_vcf_files:
	input: expand('results/tumour_sample_vcfs_octopus/{patient_id}.filtered3.vcf', patient_id=matched_somatic_patients)
	output: 'results/filtered3_joined.tsv'
	script: '../scripts/filtered4_files_joined.R'

rule collate_and_filter_vcf_archival_files:
	input: expand('results/tumour_sample_vcfs_octopus/{patient_id}.archival.filtered3.vcf', patient_id=matched_somatic_patients)
	output: 'results/archival_filtered3_joined.tsv'
	script: '../scripts/filtered4_files_joined.R'

rule collate_and_filter_vcf_relapse_files:
	input: expand('results/tumour_sample_vcfs_octopus/{patient_id}.relapse.filtered3.vcf', patient_id=matched_somatic_patients)
	output: 'results/relapse_filtered3_joined.tsv'
	script: '../scripts/filtered4_files_joined.R'

rule collate_and_filter_archival_octopus_vep_files:
	input: 
		vep_files=expand('results/tumour_sample_vcfs_octopus/{sample}.filtered.vep.vcf', sample=matched_somatic_patients),
		vcf_file='results/archival_filtered3_joined.tsv'
	output: 'results/filtered_archival_vep_calls_octopus_joined.tsv'
	script: '../scripts/collate_and_filter_vep_files_joined.R'

rule collate_and_filter_relapse_octopus_vep_files:
	input: 
		vep_files=expand('results/tumour_sample_vcfs_octopus/{sample}.filtered.vep.vcf', sample=matched_somatic_patients),
		vcf_file='results/relapse_filtered3_joined.tsv'
	output: 'results/filtered_relapse_vep_calls_octopus_joined.tsv'
	script: '../scripts/collate_and_filter_vep_files_joined.R'

rule get_archival_variants_with_two_reps:
	input: expand('results/tumour_sample_vcfs_octopus/{patient_id}.library_MAFs.vcf', patient_id=matched_somatic_patients)
	output: 'results/archival_variants_with_two_reps.tsv'
	script: '../scripts/get_num_variants_with_two_reps.R'

def get_joined_vcf (wildcards):
	test_sample_metadata = matched_somatic_metadata[(matched_somatic_metadata.fk_sample == wildcards.sample)]
	joined_vcf = expand('results/tumour_sample_vcfs_octopus/{patient_id}.filtered2.vcf', patient_id=test_sample_metadata['fk_britroc_number'])
	return(joined_vcf[0])

rule get_num_fixation_artifacts:
	input:  get_joined_vcf
	output:
		'plots/{sample}_C_to_T.png',
		'results/artifacts/{sample}.tsv'
	script: '../scripts/get_artifact_distributions.R'

rule get_num_fixation_artifact_plots:
	input: expand('results/artifacts/{sample}.tsv', sample=matched_somatic_samples)
	output: 'fixation_artifact_plot.png'
	script: '../scripts/get_artifact_plot.R' 

# maybe we should perform joint-variant calling on all archival samples belonging to the same patient?

#rule collate_and_filter_relapse_octopus_vep_files:
#	input: expand('results/tumour_sample_vcfs_octopus/{sample}.filtered.vep.vcf', sample=somatic_tp53_samples)
#	output: 'results/filtered_relapse_vep_calls_octopus.tsv'
#	script: '../scripts/collate_and_filter_vep_files.R'

#rule collate_and_filter_relapse_vep_files:
#	input: expand('results/tumour_sample_vcfs/{sample}.filtered.vep.vcf', sample=paired_relapse_samples)
#	output: 'results/filtered_relapse_vep_calls.tsv'
#	script: '../scripts/collate_and_filter_vep_files.R'

#rule collate_allele_fraction_data:
#	input: expand('results/tumour_sample_vcfs_octopus/{sample}.filtered2.vcf', sample=somatic_tp53_samples)
#	output: 'results/tp53_collated_MAFs.tsv'
#	script: '../scripts/extract_info_from_vcf.R'

#rule get_shared_variants:
#	input: expand('results/tumour_sample_vcfs/{sample}.filtered.vcf', sample=paired_somatic_samples)
#	output: 'results/shared_calls.tsv'
#	script: '../scripts/get_shared_calls.R'
