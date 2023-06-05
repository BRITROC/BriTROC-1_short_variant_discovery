library(magrittr)
library(DBI)
library(RPostgres)
library(ggplot2)
library(patchwork)

generate_somatic_oncoprint = function(somatic_variants, somatic_oncoprint_output_file, gene_set_analysed) {

	source('~/.Renviron')
	source('functions.R')

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
		dplyr::filter(!gene_symbol %in% c('RAD51L3-RFFL','CTD-2196E14.3','SHOE')) %>% dplyr::pull(gene_symbol) %>% unique()
	print(genes_with_variants)
	rm(somatic_variants_tmp)

	# remove redundant genes
	somatic_variants = somatic_variants %>% dplyr::filter(!gene_symbol %in% c('RAD51L3-RFFL','CTD-2196E14.3'))

	# remove NA entries
	somatic_variants = somatic_variants %>% dplyr::filter(!is.na(variant_type)) %>% dplyr::filter(variant_type != 'NA')

	# enforce variant type legend order and gene order
	somatic_variants$variant_type = factor(
		somatic_variants$variant_type,
		levels=c('frameshift','stop gained','splice region SNV','missense','inframe indel','no mutation')
	)
	
	somatic_variants$gene_symbol = factor(
		somatic_variants$gene_symbol,
		levels=rev(c('TP53','BRCA1','BRCA2','NF1','BRIP1','PALB2','CDK12','NRAS','RB1','KRAS','BARD1','RAD51D','RAD51C','RAD51B','PIK3CA','PTEN'))
	)

	# fill in blanks for those patient-gene combinations without a mutation
	# add in data from patients with no mutations in any gene type
	metadata = readr::read_tsv('config/all_tumour_metadata.tsv')
	somatic_variants = somatic_variants %>% dplyr::filter(patient_id %in% metadata$fk_britroc_number)	

	#no_mutation_table = tibble::tibble(
	#	patient_id=somatic_variants$patient_id %>% unique(),
	#	variant_type='no mutation'	
	#)

	no_mutation_table = tibble::tibble(
		patient_id=metadata$fk_britroc_number %>% unique(),
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

	somatic_variants$variant_type_TP53 = 
		dplyr::if_else(
			somatic_variants$gene_symbol == 'TP53',
			somatic_variants$variant_type %>% as.character(),
			'no mutation'
		)

	somatic_variants$variant_type_TP53 = somatic_variants$variant_type_TP53 %>%
		dplyr::recode(
			'frameshift'=1L,
			'missense'=1L,
			'splice region SNV'=1L,
			'stop gained'=1L,
			'inframe indel'=1L,
			'no mutation'=200L
		)

	somatic_variants$variant_type_BRCA1 = 
		dplyr::if_else(
		somatic_variants$gene_symbol == 'BRCA1',
		somatic_variants$variant_type %>% as.character(),
		'no mutation'
	)

	somatic_variants$variant_type_BRCA1 = somatic_variants$variant_type_BRCA1 %>%
		dplyr::recode(
			'frameshift'=2L,
			'missense'=2L,
			'splice region SNV'=2L,
			'stop gained'=2L,
			'inframe indel'=2L,
			'no mutation'=200L
		)

	somatic_variants$variant_type_BRCA2 = 
		dplyr::if_else(
			somatic_variants$gene_symbol == 'BRCA2',
			somatic_variants$variant_type %>% as.character(),
			'no mutation'
		)

	somatic_variants$variant_type_BRCA2 = somatic_variants$variant_type_BRCA2 %>%
		dplyr::recode(
			'frameshift'=3L,
			'missense'=3L,
			'splice region SNV'=3L,
			'stop gained'=3L,
			'inframe indel'=3L,
			'no mutation'=200L
		)

	sum_table = somatic_variants %>% dplyr::group_by(
		patient_id
	) %>% dplyr::summarise(sum=sum(variant_type_TP53, variant_type_BRCA1, variant_type_BRCA2))

	print(sum_table)
	print(somatic_variants)

	somatic_variants = 
		dplyr::left_join(
			somatic_variants,
			sum_table,
			by='patient_id'
		)

	print(somatic_variants, width=Inf)

	#quit()

	somatic_variants = somatic_variants %>% dplyr::select(-variant_type_TP53, -variant_type_BRCA1, -variant_type_BRCA2)

	somatic_variants$patient_id = 
		forcats::fct_reorder(
			somatic_variants$patient_id, 
			somatic_variants$sum, 
			.fun=min
	)

	gene_symbol_percentages = 
		somatic_variants %>%
		dplyr::filter(variant_type!='no mutation') %>%
		dplyr::group_by(gene_symbol) %>%
		dplyr::summarise(n = dplyr::n()) %>%
		dplyr::mutate(percentage = n / length(somatic_variants$patient_id %>% unique())) %>%
		.$percentage %>%
                `*`(100) %>%  
		round(digits=0) %>%
		as.character() %>%
		paste('%',sep='')

	print(gene_symbol_percentages)
	class(gene_symbol_percentages) %>% print()	

	# drop unused levels
	somatic_variants$gene_symbol = droplevels(somatic_variants$gene_symbol)
	droplevels(somatic_variants$gene_symbol) %>% levels() %>% print()
	somatic_variants$gene_symbol %>% unique() %>% print()

	p1 = ggplot(somatic_variants, aes(x=patient_id, y=as.numeric(gene_symbol))) +
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
                ) + 
		ylab('') +
                scale_y_continuous(
			breaks = 1:length(somatic_variants$gene_symbol %>% unique()),
			labels = 
				droplevels(somatic_variants$gene_symbol) %>% levels(),
				sec.axis = sec_axis(
					~.,
					breaks = 1:length(somatic_variants$gene_symbol %>% unique()),
					labels = gene_symbol_percentages
					)
		)

	patient_table = readr::read_tsv(snakemake@input['patient_table_file_path'])
	patient_table = patient_table %>% dplyr::select(britroc_number, pt_sensitivity_at_reg)
	patient_table$britroc_number = patient_table$britroc_number %>% as.factor()

	new_combined_table = dplyr::inner_join(somatic_variants, patient_table, by=c('patient_id'='britroc_number'))
	new_combined_table$pt_sensitivity_at_reg = new_combined_table$pt_sensitivity_at_reg %>% as.factor()

	new_combined_table = new_combined_table %>% dplyr::select(patient_id,pt_sensitivity_at_reg,sum) %>%
		dplyr::mutate(goo='Histological type') %>% 
		unique()

	new_combined_table$patient_id =
		forcats::fct_reorder(
			new_combined_table$patient_id,
			new_combined_table$sum,
			.fun=min
	)

	print(new_combined_table$patient_id %>% levels)
	print(somatic_variants$patient_id %>% levels %>% length())
	print(new_combined_table$patient_id %>% unique() %>% length())

	new_combined_table$pt_sensitivity_at_reg = new_combined_table$pt_sensitivity_at_reg %>%
		dplyr::recode(
			'sensitive'='Pt sensitive',
			'resistant'='Pt resistant'
		)

	p5 = ggplot(new_combined_table,
                aes(x=patient_id, y=goo, width=1.0)) +  geom_tile(aes(fill = pt_sensitivity_at_reg), colour = "white") +
                scale_fill_manual(values=c("#F8766D", "#00BFC4")) +
                theme(
                        axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        legend.title=element_blank(),
                        legend.direction='horizontal',
                        plot.margin=margin(0.001, 0.001, 0.001, 0.001, 'cm'),
                        legend.margin=margin(0, 0, 0, 0),
                        legend.box.margin=margin(-10,-10,-10,-10),
                        legend.position='top'
                ) + ylab('')

	p_final = p5 / p1 + plot_layout(heights = c(1, 16))

	ggsave(somatic_oncoprint_output_file, p_final, dev='png', width=16.0, height=9.0, scale=0.75, dpi=300)

	return()
}

generate_somatic_oncoprint(
	somatic_variants=snakemake@input[['data_for_somatic_oncoprint']],
	somatic_oncoprint_output_file=snakemake@output[['somatic_oncoprint']],
	gene_set_analysed=snakemake@params$gene_set_analysed
)
