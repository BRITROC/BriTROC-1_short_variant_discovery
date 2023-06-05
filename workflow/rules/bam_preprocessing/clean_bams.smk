# filter out suspected alignment errors in bam files by removing reads which do not map to a primer start site in the correct orientation and also by removing reads in which its mate/pair does not also map to the corresponding primer site for the same amplicon

rule symlink_bam_or_bam_indexes:
	input: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.bwamem.homo_sapiens.{bam_or_bai}'
	output: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.{bam_or_bai}'
	shell: 'ln {input} {output}'

rule create_nonoverlapping_amplicons_TP53:
	input: 'resources/union_of_tp53_amplicons.amplicons.tsv'
	output: 'resources/union_of_tp53_amplicons.amplicons.grouped.tsv',
		'resources/union_of_tp53_amplicons.targets.1.bed',
		'resources/union_of_tp53_amplicons.targets.2.bed',
		'resources/union_of_tp53_amplicons.targets.3.bed',
		'resources/union_of_tp53_amplicons.targets.4.bed',
		'resources/union_of_tp53_amplicons.targets.5.bed',
		'resources/union_of_tp53_amplicons.targets.6.bed',
		'resources/union_of_tp53_amplicons.targets.7.bed',
		'resources/union_of_tp53_amplicons.amplicons.1.bed',
		'resources/union_of_tp53_amplicons.amplicons.2.bed',
		'resources/union_of_tp53_amplicons.amplicons.3.bed',
		'resources/union_of_tp53_amplicons.amplicons.4.bed',
		'resources/union_of_tp53_amplicons.amplicons.5.bed',
		'resources/union_of_tp53_amplicons.amplicons.6.bed',
		'resources/union_of_tp53_amplicons.amplicons.7.bed'
	container: 'docker://crukcibioinformatics/ampliconseq'
	shell: 'Rscript --vanilla /opt/ampliconseq/build/bin/create_non_overlapping_amplicon_groups.R \
			--amplicons {input} \
			--reference-sequence-index {config[reference_genome_index]} \
			--output {output[0]} \
			--amplicon-bed-prefix resources/union_of_tp53_amplicons.amplicons \
			--target-bed-prefix resources/union_of_tp53_amplicons.targets'

rule create_nonoverlapping_amplicons_panel_6_28:
	input: 'resources/panel_6_28.amplicons.tsv'
	output: 'resources/panel_6_28.amplicons.grouped.tsv',
		'resources/panel_6_28.targets.1.bed',
		'resources/panel_6_28.targets.2.bed',
		'resources/panel_6_28.targets.3.bed',
		'resources/panel_6_28.targets.4.bed',
		'resources/panel_6_28.targets.5.bed',
		'resources/panel_6_28.targets.6.bed',
		'resources/panel_6_28.amplicons.1.bed',
		'resources/panel_6_28.amplicons.2.bed',
		'resources/panel_6_28.amplicons.3.bed',
		'resources/panel_6_28.amplicons.4.bed',
		'resources/panel_6_28.amplicons.5.bed',
		'resources/panel_6_28.amplicons.6.bed',
	container: 'docker://crukcibioinformatics/ampliconseq'
	shell: 'Rscript --vanilla /opt/ampliconseq/build/bin/create_non_overlapping_amplicon_groups.R \
			--amplicons {input} \
			--reference-sequence-index {config[reference_genome_index]} \
			--output {output[0]} \
			--amplicon-bed-prefix resources/panel_6_28.amplicons \
			--target-bed-prefix resources/panel_6_28.targets'

rule create_nonoverlapping_amplicons_panel_28:
	input: 'resources/panel_28.amplicons.tsv'
	output: 'resources/panel_28.amplicons.grouped.tsv',
		'resources/panel_28.targets.1.bed',
		'resources/panel_28.targets.2.bed',
		'resources/panel_28.targets.3.bed',
		'resources/panel_28.targets.4.bed',
		'resources/panel_28.targets.5.bed',
		'resources/panel_28.targets.6.bed',
		'resources/panel_28.amplicons.1.bed',
		'resources/panel_28.amplicons.2.bed',
		'resources/panel_28.amplicons.3.bed',
		'resources/panel_28.amplicons.4.bed',
		'resources/panel_28.amplicons.5.bed',
		'resources/panel_28.amplicons.6.bed',
	container: 'docker://crukcibioinformatics/ampliconseq'
	shell: 'Rscript --vanilla /opt/ampliconseq/build/bin/create_non_overlapping_amplicon_groups.R \
			--amplicons {input} \
			--reference-sequence-index {config[reference_genome_index]} \
			--output {output[0]} \
			--amplicon-bed-prefix resources/panel_28.amplicons \
			--target-bed-prefix resources/panel_28.targets'

def remove_trailing_dots(wildcards):
	mnt_bam_path = '/SLX/' + wildcards.SLX_ID + '/' + 'bam/' + wildcards.SLX_ID + '.' + wildcards.barcode + '.' + wildcards.flowcell + '.s_' + wildcards.lane + '.bam'
	return mnt_bam_path

rule clean_bams:
	input:
		bam='../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.bam',
		bai='../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.bai',
		amplicon_intervals= get_relevant_amplicon_bed_file
	output:
		bam=protected('cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.{nonoverlapping_id}.{amplicon_set_type}.bam'),
		bam_index=protected('cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.{nonoverlapping_id}.{amplicon_set_type}.bai'),
		coverage='cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.{nonoverlapping_id}.{amplicon_set_type}.coverage.txt'
	container: 'docker://crukcibioinformatics/ampliconseq'
	params:
		mnt_bam_path=remove_trailing_dots
	shell: 
		'extract-amplicon-regions \
		--id {wildcards.SLX_ID}.{wildcards.barcode}.{wildcards.flowcell}.{wildcards.lane}.{wildcards.nonoverlapping_id} \
		--unmark-duplicate-reads \
		--require-both-ends-anchored \
		--maximum-distance 1 \
		--coverage {output.coverage} \
		--input {params.mnt_bam_path} \
		--output {output.bam} \
		--amplicon-intervals {input.amplicon_intervals}'
