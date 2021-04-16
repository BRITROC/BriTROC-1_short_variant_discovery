# A snakefile within a larger workflow to generate a panel of normals for mutect2

rule link_amplicon_panels_targets_file:
	input: '/scratcha/jblab/amplicon_panels/28_JBLAB_AAprimers_dream_panel/targets.txt'
	output: 'resources/panel_28_targets.interval_list'
	shell: 'ln -s {input} {output}'

rule link_amplicon_panels_amplicon_file:
	input: '/scratcha/jblab/amplicon_panels/28_JBLAB_AAprimers_dream_panel/amplicons.txt'
	output: 'resources/panel_28_amplicons.interval_list'
	shell: 'ln -s {input} {output}'

rule mutect2_normal_only:
	input:
		reference_genome=config['reference_genome'],
		bam='../SLX/{slx}/bam/{slx}.{barcode}.{flowcell}.s_{lane}.bam',
		germline_resource='resources/gnomad.exomes.r2.1.1.fix_chr_names.sites.vcf.bgz',
		interval_file=rules.link_amplicon_panels_targets_file.output
	output: 'results/normal_sample_vcfs/{slx}_{barcode}_{flowcell}_{lane}.vcf'
	threads: 4
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk Mutect2 \
			--reference {input.reference_genome} \
			--input {input.bam} \
			--intervals {input.interval_file} \
			--germline-resource {input.germline_resource} \
			--output {output}'

rule create_list_file:
	input: expand('results/normal_sample_vcfs/{SLX}_{barcode}_{flowcell}_{lane}.vcf', zip, SLX=germline_panel_28_metadata['fk_slx'], barcode=germline_panel_28_metadata['fk_barcode'], flowcell=germline_panel_28_metadata['flowcell'], lane=germline_panel_28_metadata['lane'])
	output: 'resources/normals_for_pon_vcf.args'
	shell: 'ls {input} > {output}' 

rule create_panel_of_normals:
	input: rules.create_list_file.output
	output: 'resources/pon.vcf.gz'
	threads: 4
	shell: '/home/bioinformatics/software/gatk/gatk-4.0.11.0/gatk CreateSomaticPanelOfNormals --vcfs {input} --output {output}'

#rule normal_sample_pileup_summary:
#	input:
#		bam='../SLX/{slx}/bam/{slx}.{barcode}.{flowcell}.s_{lane}.bam',
#		biallelic_snps='somatic-b37_small_exac_common_3.fix_chr_names.vcf'
#	output: 'normal_vcfs/{slx}_{barcode}_{flowcell}_{lane}_pileups.table'
#	threads: 4
#	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk GetPileupSummaries \
#			-I {input.bam} \
#			-V {input.biallelic_snps} \
#			-L {input.biallelic_snps} \
#			-O {output}'
