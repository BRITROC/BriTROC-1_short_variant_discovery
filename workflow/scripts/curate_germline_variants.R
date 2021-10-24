library(DBI)
library(RPostgres)
library(magrittr)

curate_germline_variants = function(germline_variants_output_file) {

	source('~/.Renviron')
	source('functions.R')
	
	# make a connection to the BriTROC-1 database
	britroc_con = make_connection_to_postgres_server('britroc1', 'jblab-db.cri.camres.org', 5432)
	
	# read in germline SNVs
	germline_snvs = dbReadTable(britroc_con, 'germline_snvs') %>% 
	  dplyr::select(fk_sample,gene_symbol,chromosome,position,allele_fraction_1,allele_fraction_2,variant_type,confidence,cdna_effect,protein_effect,
	                existing_variation,clinical_signficance,thousand_genomes_allele_frequency,gnomad_allele_frequency)
	
	# read in germline indels
	germline_indels = dbReadTable(britroc_con, 'germline_indels') %>% 
	  dplyr::select(fk_sample,gene_symbol,chromosome,position,allele_fraction_1,allele_fraction_2,variant_type,confidence,cdna_effect,protein_effect,
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
	
	# generate a HGVSc column 
	germline_variants_all = germline_variants_all %>% tidyr::unite(HGVSc, gene_symbol, cdna_effect, remove=FALSE, sep=':')
	
	# filter germline variants
	germline_variants_filtered = germline_variants_all %>%
	  dplyr::filter(thousand_genomes_allele_frequency < 0.001 | is.na(thousand_genomes_allele_frequency)) %>%
	  dplyr::filter(gnomad_allele_frequency < 0.001 | is.na(gnomad_allele_frequency)) %>%
	  dplyr::filter(!grepl('benign',clinical_signficance)) %>%
	  dplyr::filter(! ( is.na(clinical_signficance) & variant_type %in% c('nonsynonymous','inframe_insertion','inframe_deletion'))) %>%
	  dplyr::filter(! ( clinical_signficance=='uncertain_significance' & variant_type %in% c('nonsynonymous','inframe_insertion', 'inframe_deletion'))) %>%
	  dplyr::filter(!variant_type %in% c('synonymous','intron','upstream_gene','downstream_gene','3_prime_UTR','5_prime_UTR'))
	
	# remove some variant as a result of extra annotation by MTBP
	final_germline_set = germline_variants_filtered %>% 
	  dplyr::filter(gene_symbol!='FANCM') %>% # extra annotation by MTBP
	  dplyr::filter(!HGVSc %in% c('BRCA1:c.4931C>G', 'BRCA1:c.5270T>C', 'BRCA1:c.5186C>A', 'RAD51B:c.453-7C>T', 'PALB2:c.1685-7T>G')) # extra annotation by MTBP
	
	# write to file
	readr::write_tsv(final_germline_set, germline_variants_output_file, append=FALSE)
	return()
	}

curate_germline_variants(snakemake@output$filtered_germline_variants)
#curate_germline_variants('results/final_germline_variants.tsv')



