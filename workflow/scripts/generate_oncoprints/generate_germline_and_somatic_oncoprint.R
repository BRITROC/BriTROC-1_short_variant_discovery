library(magrittr)
library(DBI)
library(RPostgres)
library(ComplexHeatmap)

generate_somatic_oncoprint = function(somatic_variants, germline_data, somatic_oncoprint_output_file, gene_set_analysed) {

	source('~/.Renviron')
	source('functions.R')

	# establish database connections
	britroc_con = make_connection_to_postgres_server('britroc1', 'jblab-db.cri.camres.org', 5432)

	# read in preprocessed data
	somatic_variants = readr::read_tsv(somatic_variants)

	## read in germline variant information
	germline_variants = readr::read_tsv(germline_data) %>% 
		dplyr::rename(patient_id='fk_britroc_number', Consequence='variant_type', SYMBOL='gene_symbol')
	germline_variants = germline_variants %>% dplyr::select(patient_id,Consequence,SYMBOL) %>% dplyr::mutate(type='germline')

	# only examine germline information for relevant patients
	germline_variants = germline_variants %>% dplyr::filter(patient_id %in% somatic_variants$patient_id)

	## recode germline variant information
	germline_variants$Consequence = dplyr::recode(
		germline_variants$Consequence,
		'splice_acceptor'='splice_acceptor_variant',
		'frameshift'='frameshift_variant',
		'nonsynonymous'='missense_variant',
		'stop_gained,splice_region'='stop_gained,splice_region_variant',
		'splice_region,intron'='splice_region_variant,intron_variant',
		'framehift,splice_region'='frameshift_variant'
	)

	# remap variant type categories
	germline_variants$Consequence = dplyr::recode(
	  germline_variants$Consequence,
	  'frameshift,splice_region'='frameshift',
	  'stop_gained,frameshift_variant'='frameshift',
	  'frameshift_variant'='frameshift',
	  'inframe_deletion'='inframe indel',
	  'inframe_insertion'='inframe indel',
	  'stop_gained,splice_region_variant'='stop gained',
	  'stop_gained'='stop gained',
	  'splice_acceptor'='splice region SNV',
	  'splice_region,intron'='splice region SNV',
	  'splice_region,synonymous'='splice region SNV',
	  'splice_acceptor_variant'='splice region SNV',
	  'splice_donor_variant'='splice region SNV',
	  'nonsynonymous'='missense',
	  'missense_variant'='missense',
	  "splice_region_variant,intron_variant" = 'splice region SNV',
	  "splice_region_variant,synonymous_variant"='splice region SNV',
	  'splice_region_variant'='splice region SNV',
	  'dummy'='' # make the dummy row blank
	)
	
	germline_variants = germline_variants %>% dplyr::select(-type) %>%
		dplyr::rename(gene_symbol_tumour_type=SYMBOL, variant_type=Consequence) %>%
		dplyr::select(patient_id, gene_symbol_tumour_type, variant_type)

	germline_variants_archival = germline_variants %>% dplyr::mutate(gene_symbol_tumour_type = paste(gene_symbol_tumour_type,'archival',sep='_') )
	germline_variants_relapse = germline_variants %>% dplyr::mutate(gene_symbol_tumour_type = paste(gene_symbol_tumour_type,'relapse',sep='_') )

	germline_variants = rbind(germline_variants_archival, germline_variants_relapse)
	germline_variants = germline_variants %>% dplyr::mutate(variant_type = paste(variant_type,'germline',sep=';'))

	print(somatic_variants)
	print(germline_variants)

	somatic_variants = rbind(somatic_variants,germline_variants)
	print(somatic_variants)

        # for simplicity for now - only one mutation per patient-gene group
        #somatic_variants = somatic_variants %>% tibble::as_tibble() %>%
         # dplyr::group_by(patient_id,gene_symbol_tumour_type) %>%
	 # dplyr::summarise(n=dplyr::n()) %>% dplyr::arrange(-n)

	#print(somatic_variants, n=Inf)
	#quit()
         somatic_variants = somatic_variants %>% dplyr::group_by(patient_id,gene_symbol_tumour_type) %>% 
		dplyr::filter(dplyr::row_number()==1) %>% dplyr::ungroup()

	# create a genes column
	somatic_variants_tmp = somatic_variants %>% 
		dplyr::mutate(gene=stringr::str_extract(gene_symbol_tumour_type, pattern='[A-Z0-9]+_') %>% 
		stringr::str_remove('_')) 
	
	# identify genes in which a variant has been found in at least one patient
	# is later used to remove genes in which did not have any observed variants for all tumour types
	genes_with_variants = somatic_variants_tmp %>% dplyr::pull(gene) %>% unique()
	rm(somatic_variants_tmp)

	# convert the data frame to a matrix
	somatic_variants = somatic_variants %>% 
	  tidyr::pivot_wider(names_from=patient_id, values_from=variant_type) %>% 
	  as.matrix()

	# set the row names to gene-tumour_type labels
	row.names(somatic_variants) = somatic_variants[,1]
	# remove the gene-tumour type column
	somatic_variants = somatic_variants[,-1]
	
	# replace NA values
	somatic_variants = tidyr::replace_na(somatic_variants, replace = '')
	
	# map short mutation classes to colour labels
	col = c(
	  'frameshift'= 'red',
	  'stop gained'= 'orange',
	  'splice region SNV'= 'yellow',
	  'inframe indel'= 'cyan',
	  'missense'= 'blue'
	)
	alter_fun = list(
		background = ComplexHeatmap::alter_graphic("rect", fill = "#CCCCCC"),
		frameshift = ComplexHeatmap::alter_graphic("rect", fill = col["frameshift"]),
		`stop gained` = ComplexHeatmap::alter_graphic("rect", fill = col["stop gained"]),
		`splice region SNV` = ComplexHeatmap::alter_graphic("rect", fill = col["splice region SNV"]),
		`inframe indel` = ComplexHeatmap::alter_graphic("rect", fill = col["inframe indel"]),
		missense = ComplexHeatmap::alter_graphic("rect", fill = col["missense"]),
		germline = function(x, y, w, h) grid.points(x, y, pch = 16, size = unit(0.40, "char"))
	)

	# all genes examined in this analysis
	all_genes =  tibble::tibble(names=gene_set_analysed)

	# all possible gene x tumour type combinations
	all_possible_variant_types = list(all_genes=all_genes,tumour_types=c('archival','relapse')) %>% 
		purrr::cross() %>% 
		purrr::map(purrr::lift(paste), sep='_') %>% 
		unlist()
	
	# force the all_possible_variant_types vector to have the desired order
	all_possible_variant_types = tibble::tibble(names=rep(all_genes$names,2),variant_types=all_possible_variant_types)
	print('foo')
	print(all_genes)
	print(all_possible_variant_types)
	all_genes = dplyr::inner_join(all_genes,all_possible_variant_types, by='names')
	print('goo')
	all_possible_variant_types = all_genes$variant_types

	# dummy rows to be added for gene x tumour type mutation combinations not observed
	rows_to_add_to_somatic_variants_table = setdiff(all_possible_variant_types, row.names(somatic_variants))

	dummy = matrix('', nrow=length(rows_to_add_to_somatic_variants_table), ncol=dim(somatic_variants)[2] )
	row.names(dummy) = rows_to_add_to_somatic_variants_table	

	somatic_variants = rbind(somatic_variants, dummy)

	# remove genes which have no variants
	# TODO: wrap this in a function
	somatic_variants_rows_to_keep = tibble::tibble(
		gene_symbol_tumour_type=row.names(somatic_variants),
		gene_name=stringr::str_extract(gene_symbol_tumour_type, pattern='[A-Z0-9]+_') %>% stringr::str_remove('_')
	)
	somatic_variants_rows_to_keep$keep = ifelse(somatic_variants_rows_to_keep$gene_name %in% genes_with_variants, TRUE, FALSE)

	# remove rows
	somatic_variants = somatic_variants[somatic_variants_rows_to_keep$keep,]

	# remove irrelevant items in this vector used to set row order of the oncoprint
	all_possible_variant_types = all_possible_variant_types[all_possible_variant_types %>% stringr::str_remove('_.*') %in% genes_with_variants]
	
	# generate a platinum sensitivity partition
	heatmap_patients = colnames(somatic_variants)
	heatmap_table = tibble::tibble(fk_britroc_number=as.integer(heatmap_patients))
	
	patients = dbReadTable(britroc_con, 'patients')
	
	heatmap_table = dplyr::inner_join(
	  heatmap_table, patients,
	  by=c('fk_britroc_number'='britroc_number')
	) %>% dplyr::select('fk_britroc_number','pt_sensitivity_at_reg')
	
	heatmap_table$pt_sensitivity_at_reg =
	  dplyr::recode(
	    heatmap_table$pt_sensitivity_at_reg,
	    'sensitive'='platinum sensitive',
	    'resistant'='platinum resistant'
	  )

	# set the order variants appear in the oncoprint
	variant_type_order = c('frameshift', 'stop gained', 'splice region SNV', 'inframe indel', 'missense')

	# order tables by archival TP53 variant type	
	somatic_variants = somatic_variants[, order( factor(somatic_variants[1,], levels=variant_type_order))]
	# ensure heatmap table has the same order as somatic variants
	heatmap_table = heatmap_table[match(colnames(somatic_variants) %>% as.integer, heatmap_table$fk_britroc_number),]       	

	print(somatic_variants[1,])
	#quit()

	# replace underscore in row names with a single space
	# TODO: implement a more elegant solution
	row.names(somatic_variants) = row.names(somatic_variants) %>% stringr::str_replace('_', ' ')
	all_possible_variant_types = all_possible_variant_types %>% stringr::str_replace('_', ' ')

	print(somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum resistant'])
	
	png(somatic_oncoprint_output_file, width=1200)	

        column_title_font_size = 8
	row_label_font_size = 6
	
	combined_oncoprint = ComplexHeatmap::oncoPrint(
	  mat=somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum resistant'],
	  alter_fun = alter_fun,
	  show_heatmap_legend=FALSE,
	  row_labels=rep('', length(all_possible_variant_types)),
	  column_title='Pt resistant',
          column_title_gp = gpar(fontsize = column_title_font_size, fontface='bold'),
	  pct_gp = gpar(fontsize = row_label_font_size),
	  row_names_gp = gpar(fontsize=row_label_font_size, fontface='bold'),
	  column_order = colnames(somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum resistant']),
	  row_order = all_possible_variant_types
	) +
	ComplexHeatmap::oncoPrint(
	  mat=somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum sensitive'],
	  alter_fun = alter_fun,
	  column_title='Pt sensitive',
          column_title_gp = gpar(fontsize = column_title_font_size, fontface='bold'),
	  pct_gp = gpar(fontsize = row_label_font_size),
	  row_names_gp = gpar(fontsize=row_label_font_size, fontface='bold'),
	  column_order = colnames(somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum sensitive']),
	  heatmap_legend_param = list(title='', labels_gp = gpar(fontsize = 6, fontface='bold'), ncol=5, ncol=2, by_row=TRUE),
	  row_order = all_possible_variant_types
	) 

	combined_oncoprint = ComplexHeatmap::draw(
		combined_oncoprint, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
	)
	
	dev.off()
	return()
}

#generate_somatic_oncoprint('somatic_variants_for_oncorprint.tsv', 'HRD_somatic_oncoprint.pdf', c('TP53','BRCA1','BRCA2','FANCM','BARD1','RAD51B','RAD51C','RAD51D','BRIP1','PALB2'))
generate_somatic_oncoprint(
	somatic_variants='results/data_for_somatic_oncoprint_HRD.tsv',
	germline_data=snakemake@input[['germline_data']],
	somatic_oncoprint_output_file=snakemake@output[['germline_and_somatic_oncoprint']],
	gene_set_analysed=snakemake@params$gene_set_analysed
)
