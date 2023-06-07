rule download_hs37d5_reference_genome:
	output: temp('reference_genomes/hs37d5.fa.gz')
	shell: 'wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz --output-document={output}'

rule decompress_compressed_hs37d5:
	input: rules.download_hs37d5_reference_genome.output
	output: temp('reference_genomes/hs37d5.fa')
	shell: 'gunzip {input}'

rule rename_hs37d5_chromosomes_and_scaffolds:
	input: rules.decompress_compressed_hs37d5.output
	output: temp('reference_genomes/hs37d5.with_renamed_chromosomes_and_scaffolds.fa')
	shell:  'bash workflow/scripts/reference_genomes/rename_chromosomes.sh {input} > {output}'

rule reorder_chromosomes:
	input: rules.rename_hs37d5_chromosomes_and_scaffolds.output
	output: 'reference_genomes/hs37d5.with_renamed_and_reordered_chromosomes_and_scaffolds.fa'
	conda: '../../../config/samtools.yaml'
	shell: 'bash workflow/scripts/reference_genomes/reorder_chromosomes.sh {input} {output}'

rule index_reference_genome:
	input: rules.reorder_chromosomes.output
	output: 'reference_genomes/hs37d5.with_renamed_and_reordered_chromosomes_and_scaffolds.fa.fai'
	conda: '../../../config/samtools.yaml'
	shell: 'samtools faidx {input}'

rule download_reference_genome_for_use_with_vep:
	output: 'reference_genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz'
	shell: 'wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz --output-document={output}'

rule decompress_reference_genome_for_use_with_vep:
	input: rules.download_reference_genome_for_use_with_vep.output
	output: 'reference_genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa'
	shell: 'gunzip {input}'

rule index_reference_genome_for_use_with_vep:
	input: rules.decompress_reference_genome_for_use_with_vep.output
	output: 'reference_genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa.fai'
	conda: '../../../config/samtools.yaml'
	shell: 'samtools faidx {input}'
