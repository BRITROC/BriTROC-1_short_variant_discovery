library(magrittr)
library(DBI)
library(RPostgres)
library(ggplot2)

generate_somatic_oncoprint = function(somatic_variants, somatic_oncoprint_output_file, gene_set_analysed) {

	source('~/.Renviron')
	source('functions.R')

	# establish database connections
	britroc_con = make_connection_to_postgres_server('britroc1', 'jblab-db.cri.camres.org', 5432)

	# read in preprocessed data
	somatic_variants = readr::read_tsv(somatic_variants)

	# all genes examined in this analysis
	all_genes =  tibble::tibble(names=gene_set_analysed) %>% dplyr::filter(names %in% somatic_variants$gene_symbol)

	# create a genes column
	somatic_variants_tmp = somatic_variants %>% 
		dplyr::mutate(gene=stringr::str_extract(gene_symbol, pattern='[A-Z0-9]+_') %>% 
		stringr::str_remove('_')) 
	
	# identify genes in which a variant has been found in at least one patient
	# is later used to remove genes in which did not have any observed variants for all tumour types
	genes_with_variants = somatic_variants_tmp %>% dplyr::filter(!is.na(variant_type)) %>%  
		dplyr::filter(!gene_symbol %in% c('RAD51L3-RFFL','SHOE')) %>% dplyr::pull(gene_symbol) %>% unique()
	print(genes_with_variants)
	rm(somatic_variants_tmp)

	# remove redundant genes
	somatic_variants = somatic_variants %>% dplyr::filter(!gene_symbol %in% c('RAD51L3-RFFL'))

	# remove NA entries
	somatic_variants = somatic_variants %>% dplyr::filter(!is.na(variant_type)) %>% dplyr::filter(variant_type != 'NA')

	# enforce variant type legend order and gene order
	somatic_variants$variant_type = factor(
		somatic_variants$variant_type,
		levels=c('frameshift','stop gained','splice region SNV','missense','inframe indel','no mutation')
	)
	
	somatic_variants$gene_symbol = factor(
		somatic_variants$gene_symbol,
		levels=rev(c('TP53','BRCA1','BRCA2','NF1','BRIP1','PALB2','CDK12','NRAS','RB1','KRAS','BARD1','RAD51D','RAD51B','PIK3CA','PTEN'))
	)

	# fill in blanks for those patient-gene combinations without a mutation
	no_mutation_table = tibble::tibble(
		patient_id=somatic_variants$patient_id %>% unique(),
		variant_type='no mutation'	
	)
	gene_mutation_table_tmp = tibble::tibble(
		gene_symbol=somatic_variants$gene_symbol %>% unique()
	)
	no_mutation_table = dplyr::inner_join(
		no_mutation_table,
		gene_mutation_table_tmp,
		by = character()
	)
	print(no_mutation_table)

	# remove extant mutation combinations from table
	no_mutation_table = 
		dplyr::anti_join(
			no_mutation_table,
			somatic_variants,
			by=c('patient_id','gene_symbol')
		)

	# combine both variant type tables together
	somatic_variants = rbind(
		somatic_variants,
		no_mutation_table
	)

	# convert patient_id to a factor
	somatic_variants$patient_id = somatic_variants$patient_id %>% as.factor()
	somatic_variants$gene_symbol = somatic_variants$gene_symbol %>% as.factor()

	oncoprint_data$variant_type_TP53 = 
		dplyr::if_else(
			oncoprint_data$gene_symbol == 'TP53',
			oncoprint_data$variant_type %>% as.character(),
			'no mutation'
		)

	oncoprint_data$variant_type_TP53 = oncoprint_data$variant_type_TP53 %>%
		dplyr::recode(
		'frameshift'=1L,
		'missense'=4L,
		'splice region SNV'=3L,
		'stop gained'=2L,
		'inframe indel'=5L,
		'no mutation'=200L
		)

	oncoprint_data$variant_type_BRCA1 = 
		dplyr::if_else(
		oncoprint_data$gene_symbol == 'BRCA1',
		oncoprint_data$variant_type %>% as.character(),
		'no mutation'
	)

	oncoprint_data$variant_type_BRCA1 = oncoprint_data$variant_type_BRCA1 %>%
		dplyr::recode(
			'frameshift'=5L,
			'missense'=8L,
			'splice region SNV'=7L,
			'stop gained'=6L,
			'inframe indel'=9L,
			'no mutation'=200L
		)

	oncoprint_data$variant_type_BRCA2 = 
		dplyr::if_else(
			oncoprint_data$gene_symbol == 'BRCA2',
			oncoprint_data$variant_type %>% as.character(),
			'no mutation'
		)

	oncoprint_data$variant_type_BRCA2 = oncoprint_data$variant_type_BRCA2 %>%
		dplyr::recode(
		'frameshift'=9L,
		'missense'=12L,
		'splice region SNV'=11L,
		'stop gained'=10L,
		'inframe indel'=13L,
		'no mutation'=200L
		)

	oncoprint_data$variant_type_TP53_BRCA1_BRCA2 = 
		oncoprint_data$variant_type_TP53 + 
		oncoprint_data$variant_type_BRCA1 +
		oncoprint_data$variant_type_BRCA2

	oncoprint_data$fk_britroc_number = 
		forcats::fct_reorder(
			oncoprint_data$fk_britroc_number, 
			oncoprint_data$variant_type_TP53_BRCA1_BRCA2, 
			.fun=min
	)

	p1 = ggplot(somatic_variants, aes(x=patient_id, y=gene_symbol)) +
                geom_tile(aes(fill = variant_type), colour='white', size=0.9) +
		scale_fill_manual(values=c('#604187','#FF1493','#FFD700','#00FFFF','#0FFF50','#D3D3D3')) +
                theme(
                        axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        legend.title=element_blank(),
                        axis.text.y=element_text(face="bold"),
                        plot.margin=margin(0.1, 0.5, 1.25, 0.1, 'cm'),
                        legend.margin=margin(0, 0, 0, 0),
                        legend.box.margin=margin(-10,-10,-10,-10),
                        legend.position='bottom'
                )

	ggsave(somatic_oncoprint_output_file, p1, dev='png', width=16.0, height=9.0, scale=0.75, dpi=300)

	return()
}

generate_somatic_oncoprint(
	somatic_variants=snakemake@input[['data_for_somatic_oncoprint']],
	somatic_oncoprint_output_file=snakemake@output[['somatic_oncoprint']],
	gene_set_analysed=snakemake@params$gene_set_analysed
)
