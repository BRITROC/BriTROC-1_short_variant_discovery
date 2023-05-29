rule convert_bed6_to_oct_format:
	input:  get_relevant_bed_file                                           
	output: 'resources/{analysis_type}.{nonoverlapping_id}.targets.oct'
	script: '../../../scripts/octopus_formatting/convert_bed6_to_octopus.R'

rule octopus_germline:
	input:
		interval_file=rules.convert_bed6_to_oct_format.output,
		normal_bams=cleaned_normal_bams,
		reference_genome=config['reference_genome']
	threads: 4
	output: tumour_vcf=protected('results/variant_analysis/germline/octopus_unmatched/{analysis_type}/{patient_id}.{nonoverlapping_id}.vcf')
	shell: '../octopus/bin/octopus \
			-R {input.reference_genome} \
				--allow-marked-duplicates \
				--allow-octopus-duplicates \
				--disable-downsampling \
				--annotations SB SD AF AD FRF \
				--filter-expression "QUAL < 100 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.25 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1" \
				--threads \
				--regions-file {input.interval_file} \
				-I {input.normal_bams} -o {output.tumour_vcf}'

rule filter_octopus_raw_calls_germline:
	input: 
		filtered_vcf=rules.octopus_germline.output
	output: 'results/variant_analysis/germline/{SLX_ID}.{barcodes}.{flowcell}.s_{lane}.1.panel_6_28.passed.vcf'
	shell: '/home/bioinformatics/software/bcftools/bcftools-1.10.2/bin/bcftools view -f PASS {input} | sed "/bcftools/d" > {output}'
