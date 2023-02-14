library(DBI)
library(RPostgres)
library(magrittr)

curate_germline_variants = function(germline_variants_output_file, HGVS_output_file, germline_variants_output_file_vcf) {

	source('~/.Renviron')
	source('functions.R')
	
	# make a connection to the BriTROC-1 database
	britroc_con = make_connection_to_postgres_server('britroc1', 'jblab-db.cri.camres.org', 5432)
	
	# read in germline SNVs
	germline_snvs = dbReadTable(britroc_con, 'germline_snvs') %>% 
	  dplyr::select(fk_sample,gene_symbol,chromosome,position,ref,alt,allele_fraction_1,allele_fraction_2,depth_1,depth_2,filter_1,filter_2,variant_type,confidence,cdna_effect,protein_effect,
	                existing_variation,clinical_signficance,thousand_genomes_allele_frequency,gnomad_allele_frequency)
	
	# read in germline indels
	germline_indels = dbReadTable(britroc_con, 'germline_indels') %>% 
	  dplyr::select(fk_sample,gene_symbol,chromosome,position,ref,alt,allele_fraction_1,allele_fraction_2,depth_1,depth_2,filter_1,filter_2,variant_type,confidence,cdna_effect,protein_effect,
	                existing_variation,clinical_signficance,thousand_genomes_allele_frequency,gnomad_allele_frequency)
	
	# collate all germline variants together
	germline_variants_all = rbind(germline_snvs, germline_indels)

	# filter for the relevant gene set
	germline_variants_all = germline_variants_all %>% dplyr::filter(gene_symbol %in% snakemake@params$gene_set_analysed)
	
	# read in samples table
	samples = dbReadTable(britroc_con, 'sample')
	
	# join the germline variants table and the samples table
	germline_variants_all = germline_variants_all %>% dplyr::inner_join(samples, by=c('fk_sample'='name')) %>%
	  dplyr::select(fk_britroc_number, dplyr::everything()) %>%
	  dplyr::arrange(fk_britroc_number)

	# filter germline variants
	germline_variants_filtered = germline_variants_all %>%
	  #dplyr::filter(thousand_genomes_allele_frequency < 0.001 | is.na(thousand_genomes_allele_frequency)) %>%
	  #dplyr::filter(gnomad_allele_frequency < 0.001 | is.na(gnomad_allele_frequency)) %>%
	  #dplyr::filter(!grepl('benign',clinical_signficance)) %>%
	  #dplyr::filter(! ( is.na(clinical_signficance) & variant_type %in% c('nonsynonymous','inframe_insertion','inframe_deletion'))) %>%
	  #dplyr::filter(! ( clinical_signficance=='uncertain_significance' & variant_type %in% c('nonsynonymous','inframe_insertion', 'inframe_deletion'))) %>%
	  dplyr::filter(!variant_type %in% c('upstream_gene','downstream_gene')) # 'synonymous','intron','3_prime_UTR','5_prime_UTR'

	# remove calls which are medium confidence and both calls are not no call
	germline_variants_filtered = germline_variants_filtered %>%
		dplyr::filter(!(confidence=='medium' & filter_1!='no call' & filter_2!='no call'))

	# generate a HGVSc column
	germline_variants_filtered$united_effect = dplyr::coalesce(germline_variants_filtered$protein_effect, germline_variants_filtered$cdna_effect) 
	germline_variants_filtered = germline_variants_filtered %>% tidyr::unite(HGVS_united, gene_symbol, united_effect, remove=FALSE, sep=':')
	germline_variants_filtered = germline_variants_filtered %>% tibble::as_tibble()	
	readr::write_tsv(file=HGVS_output_file,germline_variants_filtered %>% dplyr::select(HGVS_united) %>% unique(), col_names=FALSE)

	# write to file
	readr::write_tsv(germline_variants_filtered, germline_variants_output_file, append=FALSE)
	
	#convert to vcf format and write to file - needed for MTPB
	print(germline_variants_filtered)
	germline_variants_filtered_vcf = germline_variants_filtered %>% dplyr::rename(
		'#CHROM'='chromosome',
		'POS'='position',
		'REF'='ref',
		'ALT'='alt',
	)	
	germline_variants_filtered_vcf = germline_variants_filtered_vcf %>% dplyr::mutate(
		ID = '.',
		QUAL = 5000, # randomly chosen
		FILTER = 'PASS',
		INFO = 'foo'	
	)
	germline_variants_filtered_vcf = germline_variants_filtered_vcf %>% dplyr::select(`#CHROM`,POS,ID,REF,ALT,QUAL,FILTER,INFO)
	readr::write_tsv(germline_variants_filtered_vcf, germline_variants_output_file_vcf, col_names=TRUE, append=FALSE)
	
	# remove some variant as a result of extra annotation by MTBP
	#final_germline_set = germline_variants_filtered %>% 
	#  dplyr::filter(gene_symbol!='FANCM') %>% # extra annotation by MTBP
	#  dplyr::filter(!HGVS_united %in% c('BRCA1:c.4931C>G', 'BRCA1:c.5270T>C', 'BRCA1:c.5186C>A', 'RAD51B:c.453-7C>T', 'PALB2:c.1685-7T>G')) # extra annotation by MTBP
	
	return()
	}

curate_germline_variants(snakemake@output$filtered_germline_variants, snakemake@output$HGVS_output_file, snakemake@output$filtered_germline_variants_vcf)
#curate_germline_variants('results/final_germline_variants.tsv')



