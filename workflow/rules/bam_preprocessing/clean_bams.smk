#filter out suspected alignment errors in bam files by removing reads which do not map to a primer start site in the correct orientation and also by removing reads in which its mate/pair does not also map to the corresponding primer site for the same amplicon

#rule symlink_bams:
#	input: '../SLX/{SLX_ID}/bam/{SLX_ID}.{flowcell}_GRCh37_g1kp2_{lane}_{barcode}.altered_header.bam'
#	output: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.bam'
#	wildcard_constraints:
#		barcode='FLD[0-9]+',
#		flowcell="[A-Z0-9-]+",
#		lane='[1-8]'           #"((?!_bwamem).)*"
#	shell: 'ln {input} {output}'

#rule symlink_bam_indexes:
#	input: '../SLX/{SLX_ID}/bam/{SLX_ID}.{flowcell}_GRCh37_g1kp2_{lane}_{barcode}.altered_header.bai'
#	output: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.bai'
#	wildcard_constraints:
#		barcode='FLD[0-9]+',
#		flowcell="[A-Z0-9-]+",
#		lane='[1-8]'           #"((?!_bwamem).)*"
#	shell: 'ln {input} {output}'

rule symlink_bams:
	input: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.bwamem.homo_sapiens.bam'
	output: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.bam'
	wildcard_constraints:
		barcode='FLD[0-9]+',
		flowcell="[A-Z0-9-]+",
		lane='[1-8]'           #"((?!_bwamem).)*"
	shell: 'ln {input} {output}'

rule symlink_bam_indexes:
	input: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.bwamem.homo_sapiens.bai'
	output: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.bai'
	wildcard_constraints:
		barcode='FLD[0-9]+',
		flowcell="[A-Z0-9-]+",
		lane='[1-8]'           #"((?!_bwamem).)*"
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
	container: 
	shell: 'Rscript /opt/ampliconseq/build/create_non_overlapping_amplicon_groups.R \
			--amplicons {input} \
			--reference-sequence-index {config[reference_genome_index]} \
			--output {output[0]} \
			--amplicon-bed-prefix resources/union_of_tp53_amplicons.amplicons \
			--target-bed-prefix resources/union_of_tp53_amplicons.targets'

#rule create_nonoverlapping_amplicons_TP53:
#	input: 
#		amplicon_intervals='resources/union_of_tp53_amplicons.amplicons.interval_list',
#		target_intervals='resources/union_of_tp53_amplicons.targets.interval_list'
#	output: 
#		nonoverlapping_amplicon_intervals='resources/foo-union_of_tp53_amplicons.amplicons.interval_list',
#		nonoverlapping_target_intervals='resources/foo-union_of_tp53_amplicons.targets.interval_list'
#	shell: 
#		'java -Xms2G -Xmx6G \
		#		-classpath lib/ampliconseq-pipeline-0.7.2.jar:lib/commons-cli-1.4.jar:lib/commons-logging-1.2.jar:lib/htsjdk-2.19.0.jar:lib/snappy-java-1.1.4.jar \
		#		org.cruk.ampliconseq.util.CreateNonOverlappingAmpliconsAndTargets \
		#		--amplicons {input.amplicon_intervals} \
		#		--targets {input.target_intervals} \
		#		--reference {config[reference_genome]} \
		#		--output-prefix tp53.nonoverlapping \
		#		--minimum-number-of-sets 2'

rule create_nonoverlapping_amplicons2:
	input: 
		amplicon_intervals='resources/intersected_panel_6_28_amplicons.amplicons.interval_list',
		target_intervals='resources/intersected_panel_6_28_amplicons.targets.interval_list'
	output: 
		nonoverlapping_amplicon_intervals='resources/foo-panel_6_28_amplicons.amplicons.interval_list',
		nonoverlapping_target_intervals='resources/foo-panel_6_28_amplicons.targets.interval_list'
	shell: 
		'java -Xms2G -Xmx6G \
		-classpath lib/ampliconseq-pipeline-0.7.2.jar:lib/commons-cli-1.4.jar:lib/commons-logging-1.2.jar:lib/htsjdk-2.19.0.jar:lib/snappy-java-1.1.4.jar \
		org.cruk.ampliconseq.util.CreateNonOverlappingAmpliconsAndTargets \
		--amplicons {input.amplicon_intervals} \
		--targets {input.target_intervals} \
		--reference {config[reference_genome]} \
		--output-prefix panel_6_28.nonoverlapping \
		--minimum-number-of-sets 2'

rule create_nonoverlapping_amplicons_panel_28_only:
	input: 
		amplicon_intervals='resources/antijoined_panel_28_6_amplicons.amplicons.interval_list',
		target_intervals='resources/antijoined_panel_28_6_amplicons.targets.interval_list'
	output: 
		nonoverlapping_amplicon_intervals='resources/nonoverlapping_antijoined_panel_28_6_amplicons.amplicons.interval_list',
		nonoverlapping_target_intervals='resources/nonoverlapping_antijoined_panel_28_6_amplicons.targets.interval_list'
	shell: 
		'java -Xms2G -Xmx6G \
		-classpath lib/ampliconseq-pipeline-0.7.2.jar:lib/commons-cli-1.4.jar:lib/commons-logging-1.2.jar:lib/htsjdk-2.19.0.jar:lib/snappy-java-1.1.4.jar \
		org.cruk.ampliconseq.util.CreateNonOverlappingAmpliconsAndTargets \
		--amplicons {input.amplicon_intervals} \
		--targets {input.target_intervals} \
		--reference {config[reference_genome]} \
		--output-prefix nonoverlapping.antijoined_panel_28_6 \
		--minimum-number-of-sets 2'

rule create_nonoverlapping_amplicons_panel_28:
	input: 
		amplicon_intervals='/scratcha/jblab/amplicon_panels/28_JBLAB_AAprimers_dream_panel/amplicons.txt',
		target_intervals='/scratcha/jblab/amplicon_panels/28_JBLAB_AAprimers_dream_panel/targets.txt'
	output: 
		nonoverlapping_amplicon_intervals='resources/nonoverlapping.panel_28_amplicons.amplicons.interval_list',
		nonoverlapping_target_intervals='resources/nonoverlapping.panel_28_amplicons.targets.interval_list'
	shell: 
		'java -Xms2G -Xmx6G \
		-classpath lib/ampliconseq-pipeline-0.7.2.jar:lib/commons-cli-1.4.jar:lib/commons-logging-1.2.jar:lib/htsjdk-2.19.0.jar:lib/snappy-java-1.1.4.jar \
		org.cruk.ampliconseq.util.CreateNonOverlappingAmpliconsAndTargets \
		--amplicons {input.amplicon_intervals} \
		--targets {input.target_intervals} \
		--reference {config[reference_genome]} \
		--output-prefix nonoverlapping.panel_28 \
		--minimum-number-of-sets 2'

def remove_trailing_dots(wildcards):
	mnt_bam_path = '/SLX/' + wildcards.SLX_ID + '/' + 'bam/' + wildcards.SLX_ID + '.' + wildcards.barcode + '.' + wildcards.flowcell + '.s_' + wildcards.lane + '.bam'
	return mnt_bam_path

rule clean_bams:
	input:
		bam='../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.bam',
		bai='../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.s_{lane}.bai',
		amplicon_intervals='resources/union_of_tp53_amplicons.amplicons.{nonoverlapping_id}.bed'
	output:
		bam=protected('cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.union_of_tp53.bam'),
		bam_index=protected('cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.union_of_tp53.bai'),
		coverage='cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.union_of_tp53.coverage.txt'
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
