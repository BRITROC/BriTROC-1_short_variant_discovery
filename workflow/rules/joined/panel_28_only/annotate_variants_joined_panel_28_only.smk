matched_somatic_metadata_panel_28_only = pandas.read_table("config/somatic_metadata_panel_28_only.tsv").set_index("fk_britroc_number", drop=False)
matched_somatic_patients_panel_28_only = matched_somatic_metadata_panel_28_only.index.unique().tolist()

rule vep_octopus:
	input: rules.concat_vcfs.output
	output: 'results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.filtered.vep.vcf'
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

rule ensure_tech_rep_genotypes_match:
	input: rules.concat_vcfs.output
	output: 
		tumour_samples_union='results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.filtered3.vcf',
		archival_samples='results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.archival.filtered3.vcf',
		relapse_samples='results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.relapse.filtered3.vcf',
		library_MAFs='results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.library_MAFs.vcf',
		library_depths='results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.library_depths.vcf',
		sample_genotypes='results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.sample_genotypes.vcf'
	script: '../scripts/annotate_variants_joined/view_square_vcfs.R'

rule collate_and_filter_vcf_files:
	input: expand('results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.filtered3.vcf', patient_id=matched_somatic_patients_panel_28_only)
	output: 'results/panel_28_only/filtered3_joined.tsv'
	script: '../scripts/annotate_variants_joined/filtered4_files_joined.R'

rule collate_and_filter_vcf_archival_files:
	input: expand('results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.archival.filtered3.vcf', patient_id=matched_somatic_patients_panel_28_only)
	output: 'results/panel_28_only/archival_filtered3_joined.tsv'
	script: '../scripts/annotate_variants_joined/filtered4_files_joined.R'

rule collate_and_filter_vcf_relapse_files:
	input: expand('results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.relapse.filtered3.vcf', patient_id=matched_somatic_patients_panel_28_only)
	output: 'results/panel_28_only/relapse_filtered3_joined.tsv'
	script: '../scripts/annotate_variants_joined/filtered4_files_joined.R'

rule collate_and_filter_archival_octopus_vep_files:
	input: 
		vep_files=expand('results/tumour_sample_vcfs_octopus/panel_28_only/{sample}.filtered.vep.vcf', sample=matched_somatic_patients_panel_28_only),
		vcf_file=rules.collate_and_filter_vcf_archival_files.output
	output: 'results/panel_28_only/filtered_archival_vep_calls_octopus_joined.tsv'
	script: '../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'

rule collate_archival_octopus_vep_files_without_filtering:
	input: 
		vep_files=expand('results/tumour_sample_vcfs_octopus/panel_28_only/{sample}.filtered.vep.vcf', sample=matched_somatic_patients_panel_28_only),
		vcf_file=rules.collate_and_filter_vcf_archival_files.output
	output: 'results/panel_28_only/collated_archival_vep_calls_octopus_joined.tsv'
	script: '../scripts/annotate_variants_joined/collate_files_without_filtering.R'

rule collate_and_filter_relapse_octopus_vep_files:
	input: 
		vep_files=expand('results/tumour_sample_vcfs_octopus/panel_28_only/{sample}.filtered.vep.vcf', sample=matched_somatic_patients_panel_28_only),
		vcf_file=rules.collate_and_filter_vcf_relapse_files.output
	output: 'results/panel_28_only/filtered_relapse_vep_calls_octopus_joined.tsv'
	script: '../scripts/annotate_variants_joined/collate_and_filter_vep_files_joined.R'

rule collate_relapse_octopus_vep_files_without_filtering:
	input: 
		vep_files=expand('results/tumour_sample_vcfs_octopus/panel_28_only/{sample}.filtered.vep.vcf', sample=matched_somatic_patients_panel_28_only),
		vcf_file=rules.collate_and_filter_vcf_relapse_files.output
	output: 'results/panel_28_only/collated_relapse_vep_calls_octopus_joined.tsv'
	script: '../scripts/annotate_variants_joined/collate_files_without_filtering.R'

rule get_archival_variants_with_two_reps:
	input: expand('results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.library_MAFs.vcf', patient_id=matched_somatic_patients_panel_28_only)
	output: 'results/panel_28_only/archival_variants_with_two_reps.tsv'
	script: '../scripts/annotate_variants_joined/get_num_variants_with_two_reps.R'

def get_joined_vcf (wildcards):
	test_sample_metadata = matched_somatic_metadata_panel_28_only[(matched_somatic_metadata_panel_28_only.fk_sample == wildcards.sample)]
	joined_vcf = expand('results/tumour_sample_vcfs_octopus/panel_28_only/{patient_id}.filtered2.vcf', patient_id=test_sample_metadata_panel_28_only['fk_britroc_number'])
	return(joined_vcf[0])

#rule get_num_fixation_artifacts:
#	input:  get_joined_vcf
#	output:
#		'plots/{sample}_C_to_T.png',
#		'results/artifacts/{sample}.tsv'
#	script: '../scripts/get_artifact_distributions.R'

#rule get_num_fixation_artifact_plots:
#	input: expand('results/artifacts/{sample}.tsv', sample=matched_somatic_samples)
#	output: 'fixation_artifact_plot.png'
#	script: '../scripts/get_artifact_plot.R' 
