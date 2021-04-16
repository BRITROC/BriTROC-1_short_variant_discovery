def get_contamination_tables(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	metadata_filt = paired_somatic_metadata[(paired_somatic_metadata.fk_sample==wildcards.sample)].head(n=2)
	SLX = metadata_filt.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	barcodes = metadata_filt.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()
	lane = metadata_filt.set_index('lane', drop=False)
	lane = lane.index.unique().tolist()
	flowcell = metadata_filt.set_index('fk_run', drop=False)
	flowcell = flowcell.index.unique().tolist()
	flowcell = flowcell[0].split('_')

	contamination_table_1 = 'results/tumour_sample_vcfs/{}_{}_{}_{}_contamination.table'.format(SLX[0], barcodes[0], flowcell[-1], lane[0])
	contamination_table_2 = 'results/tumour_sample_vcfs/{}_{}_{}_{}_contamination.table'.format(SLX[0], barcodes[1], flowcell[-1], lane[0])

	return([contamination_table_1, contamination_table_2])

def get_segmentation_tables(wildcards):
	# some samples have been sequenced multiple times which is a variable we will have to factor in later
	metadata_filt = paired_somatic_metadata[(paired_somatic_metadata.fk_sample==wildcards.sample)].head(n=2)
	SLX = metadata_filt.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	barcodes = metadata_filt.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()
	lane = metadata_filt.set_index('lane', drop=False)
	lane = lane.index.unique().tolist()
	flowcell = metadata_filt.set_index('fk_run', drop=False)
	flowcell = flowcell.index.unique().tolist()
	flowcell = flowcell[0].split('_')

	segmentation_table_1 = 'results/tumour_sample_vcfs/{}_{}_{}_{}_segmentation.table'.format(SLX[0], barcodes[0], flowcell[-1], lane[0])
	segmentation_table_2 = 'results/tumour_sample_vcfs/{}_{}_{}_{}_segmentation.table'.format(SLX[0], barcodes[1], flowcell[-1], lane[0])

	return([segmentation_table_1, segmentation_table_2])

rule LearnReadOrientationModel:
	input:  rules.mutect2.output.f1r2   #'tumour_vcfs/{sample}_f1r2.tar.gz'
	output: 'results/tumour_sample_vcfs/{sample}_artifiact_priors.tar.gz'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk LearnReadOrientationModel \
			--input {input} \
			--output {output}'

rule get_tumour_pileup_summary:
	input:
		bam='../SLX/{slx}/bam/{slx}.{barcode}.{flowcell}.s_{lane}.bam',
		biallelic_snps='resources/somatic-b37_small_exac_common_3.fix_chr_names.vcf'
	output: 'results/tumour_sample_vcfs/{slx}_{barcode}_{flowcell}_{lane}_pileups.table'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk GetPileupSummaries \
			-I {input.bam} \
			-V {input.biallelic_snps} \
			-L {input.biallelic_snps} \
			-O {output}'

rule get_contamination_table:
	input: rules.get_tumour_pileup_summary.output,
	output: 
		contamination_table='results/tumour_sample_vcfs/{slx}_{barcode}_{flowcell}_{lane}_contamination.table',
		tumour_segmentation='results/tumour_sample_vcfs/{slx}_{barcode}_{flowcell}_{lane}_segmentation.table'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk CalculateContamination \
			-I {input} \
			--tumor-segmentation {output.tumour_segmentation} \
			-O {output.contamination_table}'


rule mark_mutect_calls_for_filtering:
	input:
		reference=config['reference_genome'],
		mutect2_vcf=rules.mutect2.output.tumour_vcf,
		f1r2=rules.LearnReadOrientationModel.output,
		contamination_tables=get_contamination_tables,
		segmentation_tables=get_segmentation_tables
	output: 'results/tumour_sample_vcfs/{sample}.filtered.vcf'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk FilterMutectCalls \
			-R {input.reference} \
			-V {input.mutect2_vcf} \
			--orientation-bias-artifact-priors {input.f1r2} \
			--contamination-table {input.contamination_tables[0]} \
			--contamination-table {input.contamination_tables[1]} \
			--tumor-segmentation {input.segmentation_tables[0]} \
			--tumor-segmentation {input.segmentation_tables[1]} \
			--min-allele-fraction 0.00 \
			-O {output}'

rule filter_variants:
	input: 
		filtered_vcf=rules.mark_mutect_calls_for_filtering.output,
		reference=config['reference_genome']
	output: 'results/tumour_sample_vcfs/{sample}.filtered2.vcf'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk SelectVariants \
			-R {input.reference} \
			-V {input.filtered_vcf} \
			--exclude-filtered \
			-O {output}'

rule filter_vardict_variants:
	input: 
		filtered_vcf=rules.vardict.output,
		reference=config['reference_genome']
	output: 'results/tumour_sample_vcfs_vardict/{sample}.filtered.vcf'
	shell: '/home/bioinformatics/software/gatk/gatk-4.1.8.0/gatk SelectVariants \
			-R {input.reference} \
			-V {input.filtered_vcf} \
			--exclude-filtered \
			-O {output}'

# note: Can't use SelectVariant with octopus output as they are not compatible

rule filter_octopus_raw_calls:
	input: 
		filtered_vcf='results/tumour_sample_vcfs_octopus/{sample}.vcf',
	output: 'results/tumour_sample_vcfs_octopus/{sample}.filtered.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

rule filter_octopus_hard_filtering:
	input: 
		filtered_vcf='results/tumour_sample_vcfs_octopus/{sample}.filtered2.vcf',
	output: 'results/tumour_sample_vcfs_octopus/{sample}.filtered3.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

rule filter_octopus_random_forest:
	input: 
		filtered_vcf='results/tumour_sample_vcfs_octopus/{sample}.filtered4.vcf',
	output: 'results/tumour_sample_vcfs_octopus/{sample}.filtered5.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'

# filters VCF and preserves record in which all samples (except the normal) have the same predicted genotype
rule ensure_matching_genotypes:
	input: 'results/tumour_sample_vcfs_octopus/{sample}.filtered6.vcf'
	output: 'results/tumour_sample_vcfs_octopus/{sample}.filtered7.vcf'
	shell: 'bash workflow/scripts/match_genotypes.sh {input} > {output}'
