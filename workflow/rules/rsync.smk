# A snakemake rules file to transfer BAM files between servers where thry can be viewed by applications such as IGV more easily.
# For ease of interpretation, BAM files have also been renamed - incorporating important sample metadata information in the sample name

def get_bam_and_bam_index(wildcards):
	#patient_metadata = all_somatic_metadata[(all_somatic_metadata.fk_sample == wildcards.sample_id)] 
	#patient_metadata = patient_metadata[(patient_metadata.fk_barcode == wildcards.barcode)]

	patient_metadata = germline_metadata[(germline_metadata.fk_sample == wildcards.sample_id)] 
	patient_metadata = patient_metadata[(patient_metadata.fk_barcode == wildcards.barcode)]

	SLX = patient_metadata.set_index('fk_slx', drop=False)
	SLX = SLX.index.unique().tolist()
	barcodes = patient_metadata.set_index('fk_barcode', drop=False)
	barcodes = barcodes.index.unique().tolist()
	lane = patient_metadata.set_index('lane', drop=False)
	lane = lane.index.unique().tolist()
	flowcell = patient_metadata.set_index('flowcell', drop=False)
	flowcell = flowcell.index.unique().tolist()
	flowcell = flowcell[0].split('_')

	bam_file = '../SLX/{}/bam/{}.{}.{}.s_{}.bam'.format(SLX[0],  SLX[0], barcodes[0], flowcell[0], lane[0])
	bam_index_file = '../SLX/{}/bam/{}.{}.{}.s_{}.bai'.format(SLX[0],  SLX[0], barcodes[0], flowcell[0], lane[0])

	return([bam_file, bam_index_file])

rule rsync_bams:
	input: get_bam_and_bam_index
	output: 'results/rsync/patient-{fk_britroc_number}_{sample_type}_{sample_id}_{barcode}'
	wildcard_constraints:
		sample_type='(archival|relapse|germline)'
	shell: 'rsync {input[0]} jblab-srv001:/scratch/britroc1_manuscript2/patient_bams/{wildcards.fk_britroc_number}/britroc1-patient{wildcards.fk_britroc_number}_{wildcards.sample_type}_{wildcards.sample_id}_{wildcards.barcode}.bam && rsync {input[1]} jblab-srv001:/scratch/britroc1_manuscript2/patient_bams/{wildcards.fk_britroc_number}/britroc1-patient{wildcards.fk_britroc_number}_{wildcards.sample_type}_{wildcards.sample_id}_{wildcards.barcode}.bai && touch {output}'
