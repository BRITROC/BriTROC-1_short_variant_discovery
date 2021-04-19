# filter out suspected alignment errors in bam files by removing reads which do not map to a primer start site in the correct orientation and also by removing reads in which its mate/pair does not also map to the corresponding primer site for the same amplicon

rule clean_bams:
	input:
		bam='../SLX/{SLX_ID}/bam/{SLX_ID}.{barcode}.{flowcell}.{lane}.bam',
		amplicon_intervals='resources/intersected_panel_6_28_amplicons.amplicons.interval_list'
	output:
		bam='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.amplicon.bam',
		bam_index='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.amplicon.bai',
		coverage='../SLX/{SLX_ID}/bam/cleaned_bams/{SLX_ID}.{barcode}.{flowcell}.{lane}.amplicon_coverage.txt'
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

# Apply the same logic to merged bam files which have slightly different file patterns		
rule clean_sample_bams:
	input:
		bam='../SLX/{SLX_ID}/samplebam/{SLX_ID}.{barcode}.bam',
		amplicon_intervals='resources/intersected_panel_6_28_amplicons.amplicons.interval_list'
	output:
		bam= '../SLX/{SLX_ID}/samplebam/cleaned_samplebams/{SLX_ID}.{barcode}.amplicon.bam',
		bam_index= '../SLX/{SLX_ID}/samplebam/cleaned_samplebams/{SLX_ID}.{barcode}.amplicon.bai',
		coverage= '../SLX/{SLX_ID}/samplebam/cleaned_samplebams/{SLX_ID}.{barcode}.amplicon_coverage.txt'
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
