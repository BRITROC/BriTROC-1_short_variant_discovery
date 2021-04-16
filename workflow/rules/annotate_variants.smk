rule filter_variants2:
	input: 
		filtered_vcf=rules.filter_variants.output
	output: 'results/tumour_sample_vcfs/{sample}.filtered3.vcf'
	shell: 'sed s/^chr//g {input.filtered_vcf} > {output}'

rule vep:
	input:
		filtered_vcf=rules.filter_variants.output
	output: 'results/tumour_sample_vcfs/{sample}.filtered.vep.vcf'
	shell: '/home/bioinformatics/software/ensembl/ensembl-vep-release-102.0/vep \
			-i {input} \
			-o {output} \
			--cache \
			--offline \
			--format vcf \
			--dir /mnt/scratchb/bioinformatics/reference_data/ensembl/vep \
			--force_overwrite \
			--check_existing \
			--everything \
			--no_escape \
			-a GRCh37 \
			--port 3337'

rule vep_vardict:
	input:
		filtered_vcf=rules.filter_vardict_variants.output
	output: 'results/tumour_sample_vcfs_vardict/{sample}.filtered.vep.vcf'
	shell: '/home/bioinformatics/software/ensembl/ensembl-vep-release-102.0/vep \
			-i {input} \
			-o {output} \
			--cache \
			--offline \
			--format vcf \
			--dir /mnt/scratchb/bioinformatics/reference_data/ensembl/vep \
			--force_overwrite \
			--check_existing \
			--everything \
			--no_escape \
			-a GRCh37 \
			--port 3337'

rule vep_octopus:
	input:
		'results/tumour_sample_vcfs_octopus/{sample}.filtered7.vcf'
	output: 'results/tumour_sample_vcfs_octopus/{sample}.filtered.vep.vcf'
	shell: 'ensembl-vep/vep \
			-i {input} \
			-o {output} \
			--cache \
			--offline \
			--format vcf \
			--dir /Users/bradle02/.vep/ \
			--force_overwrite \
			--check_existing \
			--everything \
			--no_escape \
			-a GRCh37 \
			--port 3337'

rule link_vep_filtering_script:
	input: '../strelka/strelka_binary/collate_and_filter_vep_files.R'
	output: 'workflow/scripts/collate_and_filter_vep_files.R'
	shell: 'ln -s {input} {output}'

rule collate_and_filter_archival_vep_files:
	input: expand('results/tumour_sample_vcfs/{sample}.filtered.vep.vcf', sample=paired_archival_samples)
	output: 'results/filtered_archival_vep_calls.tsv'
	script: '../scripts/collate_and_filter_vep_files.R'

rule collate_and_filter_archival_vardict_vep_files:
	input: expand('results/tumour_sample_vcfs_vardict/{sample}.filtered.vep.vcf', sample=paired_archival_samples)
	output: 'results/filtered_archival_vep_calls_vardict.tsv'
	script: '../scripts/collate_and_filter_vep_files.R'

rule collate_and_filter_archival_octopus_vep_files:
	input: expand('results/tumour_sample_vcfs_octopus/{sample}.filtered.vep.vcf', sample=paired_archival_samples)
	output: 'results/filtered_archival_vep_calls_octopus.tsv'
	script: '../scripts/collate_and_filter_vep_files.R'

rule collate_and_filter_relapse_octopus_vep_files:
	input: expand('results/tumour_sample_vcfs_octopus/{sample}.filtered.vep.vcf', sample=paired_relapse_samples)
	output: 'results/filtered_relapse_vep_calls_octopus.tsv'
	script: '../scripts/collate_and_filter_vep_files.R'

rule collate_and_filter_relapse_vep_files:
	input: expand('results/tumour_sample_vcfs/{sample}.filtered.vep.vcf', sample=paired_relapse_samples)
	output: 'results/filtered_relapse_vep_calls.tsv'
	script: '../scripts/collate_and_filter_vep_files.R'

rule get_shared_variants:
	input: expand('results/tumour_sample_vcfs/{sample}.filtered.vcf', sample=paired_somatic_samples)
	output: 'results/shared_calls.tsv'
	script: '../scripts/get_shared_calls.R'
