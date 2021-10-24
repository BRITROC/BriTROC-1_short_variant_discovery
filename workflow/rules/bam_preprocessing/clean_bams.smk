# filter out suspected alignment errors in bam files by removing reads which do not map to a primer start site in the correct orientation and also by removing reads in which its mate/pair does not also map to the corresponding primer site for the same amplicon

rule symlink_bams:
	input: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.{lane}.bwamem.homo_sapiens.bam'
	output: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.{lane}.bam'
	wildcard_constraints:
		barcode='FLD[0-9]+',
		flowcell="[A-Z0-9-]+",
		lane='s_[1-8]'           #"((?!_bwamem).)*"
	shell: 'ln {input} {output}'

rule symlink_bam_indexes:
	input: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.{lane}.bwamem.homo_sapiens.bai'
	output: '../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.{lane}.bai'
	wildcard_constraints:
		barcode='FLD[0-9]+',
		flowcell="[A-Z0-9-]+",
		lane='s_[1-8]'           #"((?!_bwamem).)*"
	shell: 'ln {input} {output}'

rule create_nonoverlapping_amplicons_TP53:
	input: 
		amplicon_intervals='resources/union_of_tp53_amplicons.amplicons.interval_list',
		target_intervals='resources/union_of_tp53_amplicons.targets.interval_list'
	output: 
		nonoverlapping_amplicon_intervals='resources/foo-union_of_tp53_amplicons.amplicons.interval_list',
		nonoverlapping_target_intervals='resources/foo-union_of_tp53_amplicons.targets.interval_list'
	shell: 
		'java -Xms2G -Xmx6G \
		-classpath lib/ampliconseq-pipeline-0.7.2.jar:lib/commons-cli-1.4.jar:lib/commons-logging-1.2.jar:lib/htsjdk-2.19.0.jar:lib/snappy-java-1.1.4.jar \
		org.cruk.ampliconseq.util.CreateNonOverlappingAmpliconsAndTargets \
		--amplicons {input.amplicon_intervals} \
		--targets {input.target_intervals} \
		--reference {config[reference_genome]} \
		--output-prefix tp53.nonoverlapping \
		--minimum-number-of-sets 2'

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

rule clean_bams:
	input:
		bam='../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.{lane}.bam',
		bai='../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.{lane}.bai',
		amplicon_intervals='resources/tp53.nonoverlapping.amplicons.{nonoverlapping_id}.bed'
	output:
		bam='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.union_of_tp53.bam',
		bam_index='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.union_of_tp53.bai',
		coverage='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.union_of_tp53.coverage.txt'
	wildcard_constraints:
		barcode='FLD[0-9]+'
	shell: 
		'java -Xms2G -Xmx6G \
		-classpath lib/ampliconseq-pipeline-0.7.2.jar:lib/commons-cli-1.4.jar:lib/commons-logging-1.2.jar:lib/htsjdk-2.19.0.jar:lib/snappy-java-1.1.4.jar \
		org.cruk.ampliconseq.util.ExtractAmpliconRegions \
		--unmark-duplicate-reads \
		--require-both-ends-anchored \
		--maximum-distance 1 \
		--amplicon-coverage {output.coverage} \
		--bam {input.bam} \
		--amplicon-bam {output.bam} \
		--amplicons {input.amplicon_intervals}'

rule clean_bams2:
	input:
		bam='../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.{lane}.bam',
		bai='../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.{lane}.bai',
		amplicon_intervals='resources/panel_6_28.nonoverlapping.amplicons.{nonoverlapping_id}.bed'
	output:
		bam='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.panel_6_28.bam',
		bam_index='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.panel_6_28.bai',
		coverage='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.panel_6_28.coverage.txt'
	wildcard_constraints:
		barcode='FLD[0-9]+'
	shell: 
		'java -Xms2G -Xmx6G \
		-classpath lib/ampliconseq-pipeline-0.7.2.jar:lib/commons-cli-1.4.jar:lib/commons-logging-1.2.jar:lib/htsjdk-2.19.0.jar:lib/snappy-java-1.1.4.jar \
		org.cruk.ampliconseq.util.ExtractAmpliconRegions \
		--unmark-duplicate-reads \
		--require-both-ends-anchored \
		--maximum-distance 1 \
		--amplicon-coverage {output.coverage} \
		--bam {input.bam} \
		--amplicon-bam {output.bam} \
		--amplicons {input.amplicon_intervals}'

rule clean_bams_panel_28_only:
	input:
		bam='../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.{lane}.bam',
		bai='../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.{lane}.bai',
		amplicon_intervals='resources/nonoverlapping.antijoined_panel_28_6.amplicons.{nonoverlapping_id}.bed'
	output:
		bam='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.antijoined_panel_28_6.bam',
		bam_index='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.antijoined_panel_28_6.bai',
		coverage='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.{nonoverlapping_id}.antijoined_panel_28_6.coverage.txt'
	wildcard_constraints:
		barcode='FLD[0-9]+'
	shell: 
		'java -Xms2G -Xmx6G \
		-classpath lib/ampliconseq-pipeline-0.7.2.jar:lib/commons-cli-1.4.jar:lib/commons-logging-1.2.jar:lib/htsjdk-2.19.0.jar:lib/snappy-java-1.1.4.jar \
		org.cruk.ampliconseq.util.ExtractAmpliconRegions \
		--unmark-duplicate-reads \
		--require-both-ends-anchored \
		--maximum-distance 1 \
		--amplicon-coverage {output.coverage} \
		--bam {input.bam} \
		--amplicon-bam {output.bam} \
		--amplicons {input.amplicon_intervals}'
