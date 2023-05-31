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

rule collate_archival_octopus_vep_files_without_filtering:
	input: 
		vep_files=expand('results/tumour_sample_vcfs_octopus/{sample}.filtered.vep.vcf', sample=matched_somatic_patients),
		vcf_file='results/archival_filtered3_joined.tsv'
	output: 'results/collated_archival_vep_calls_octopus_joined.tsv'
	script: '../scripts/annotate_variants_joined/collate_files_without_filtering.R'

rule collate_relapse_octopus_vep_files_without_filtering:
	input: 
		vep_files=expand('results/tumour_sample_vcfs_octopus/{sample}.filtered.vep.vcf', sample=matched_somatic_patients),
		vcf_file='results/relapse_filtered3_joined.tsv'
	output: 'results/collated_relapse_vep_calls_octopus_joined.tsv'
	script: '../scripts/annotate_variants_joined/collate_files_without_filtering.R'

rule get_archival_variants_with_two_reps:
	input: expand('results/tumour_sample_vcfs_octopus/{patient_id}.library_MAFs.vcf', patient_id=matched_somatic_patients)
	output: 'results/archival_variants_with_two_reps.tsv'
	script: '../scripts/annotate_variants_joined/get_num_variants_with_two_reps.R'

rule collate_and_filter_vcf_files:
	input: expand('results/tumour_sample_vcfs_octopus/{analysis_type}/{patient_id}.filtered3.vcf', patient_id=matched_somatic_patients)
	output: 'results/filtered3_joined.tsv'
	script: '../scripts/annotate_variants_joined/filtered4_files_joined.R'


