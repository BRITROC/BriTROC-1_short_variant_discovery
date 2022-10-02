library(magrittr)
library(DBI)
library(RPostgres)
library(ComplexHeatmap)

generate_somatic_oncoprint = function(somatic_variants, somatic_oncoprint_output_file, gene_set_analysed) {

	source('~/.Renviron')
	source('functions.R')

	# establish database connections
	britroc_con = make_connection_to_postgres_server('britroc1', 'jblab-db.cri.camres.org', 5432)

	# read in preprocessed data
	somatic_variants = readr::read_tsv(somatic_variants)

	# remove bad gene
	somatic_variants = somatic_variants %>% dplyr::filter(!grepl('RAD51L3-RFFL',gene_symbol_tumour_type))

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
	all_genes = dplyr::inner_join(all_genes,all_possible_variant_types, by='names')
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
	somatic_variants = somatic_variants[, order( factor(somatic_variants[3,], levels=variant_type_order))]
	# ensure heatmap table has the same order as somatic variants
	heatmap_table = heatmap_table[match(colnames(somatic_variants) %>% as.integer, heatmap_table$fk_britroc_number),]       	

	# replace underscore in row names with a single space
	# TODO: implement a more elegant solution
	row.names(somatic_variants) = row.names(somatic_variants) %>% stringr::str_replace('_', ' ')
	all_possible_variant_types = all_possible_variant_types %>% stringr::str_replace('_', ' ')

	print(somatic_variants)

	#TODO : remove this quick fix/fudge
	#somatic_variants = somatic_variants[-9,]

	#print(somatic_variants)
	print(all_possible_variant_types)

	print(somatic_variants)
	readr::write_tsv(somatic_variants %>% as.data.frame(), file="whole_cohort_intercalated_matrix.txt", col_names=TRUE)
	
	png(somatic_oncoprint_output_file, width=1200)	

        column_title_font_size = 8
	row_label_font_size = 6

	somatic_variants = somatic_variants[!rownames(somatic_variants) %in% 'SHOE foo', ]  # ! is logical negation

	somatic_variants %>% dim() %>% print()

	all_possible_variant_types = c(
	'TP53 archival','TP53 relapse',
	'BRCA1 archival','BRCA1 relapse',
	'BRCA2 archival','BRCA2 relapse',
	"RAD51D archival" ,"RAD51D relapse",
	"BRIP1 archival","BRIP1 relapse",
	"PALB2 archival","PALB2 relapse",
	"BARD1 archival","BARD1 relapse",
	"CDK12 archival","CDK12 relapse",
	"KRAS archival","KRAS relapse",
	"PIK3CA archival","PIK3CA relapse",
	"NF1 archival","NF1 relapse",
	"RB1 archival","RB1 relapse",
	"NRAS archival","NRAS relapse"
	)


	#print(all_possible_variant_types)
	rownames_somatic_variants = row.names(somatic_variants)
	print(rownames_somatic_variants)
	somatic_variants = somatic_variants %>% 
		tibble::as_tibble() %>% 
		dplyr::mutate(variant_type=rownames_somatic_variants) %>% 
		dplyr::select(variant_type, dplyr::everything())

	#print(somatic_variants)

	transform_data = function(gene_name) {

	# subset rows	
	somatic_variants_one_gene = somatic_variants %>% dplyr::filter(grepl(gene_name,variant_type))

	# coalesce columns
	coalesce_by_column <- function(df) {
  		return(coalesce(df[1], df[2]))
	}

	somatic_variants_one_gene_coalesced = dplyr::if_else(somatic_variants_one_gene[1,]=='' & somatic_variants_one_gene[2,]=='', '', as.character(NA))
	somatic_variants_one_gene_coalesced = dplyr::if_else(somatic_variants_one_gene[1,]!='' & somatic_variants_one_gene[2,]!='', 'shared mutation', somatic_variants_one_gene_coalesced)
	somatic_variants_one_gene_coalesced = dplyr::if_else(somatic_variants_one_gene[1,]!='' & somatic_variants_one_gene[2,]=='', 'archival only mutation', somatic_variants_one_gene_coalesced)
	somatic_variants_one_gene_coalesced = dplyr::if_else(somatic_variants_one_gene[1,]=='' & somatic_variants_one_gene[2,]!='', 'relapse only mutation', somatic_variants_one_gene_coalesced)

	#print(somatic_variants_one_gene)

	#print(somatic_variants_one_gene_coalesced)	

	somatic_variants_one_gene_coalesced[1] = gene_name
	somatic_variants_one_gene[1,] = as.list(somatic_variants_one_gene_coalesced)

	return(somatic_variants_one_gene[1,])

	}

	gene_list = all_possible_variant_types %>% stringr::str_remove_all(' archival') %>% stringr::str_remove_all(' relapse') %>% unique()
	new_somatic_variants = purrr::map_dfr(.x=gene_list, .f=transform_data)
	new_somatic_variants = new_somatic_variants %>% dplyr::select(-variant_type)
		
	print(new_somatic_variants)
	readr::write_tsv(new_somatic_variants, 'tmp.tsv')

	new_somatic_variants = as.matrix(new_somatic_variants)
	row.names(new_somatic_variants) = gene_list

	print(new_somatic_variants)

# [1] "BRCA1 archival"  "BRCA1 relapse"   "BRCA2 archival"  "BRCA2 relapse"
# [5] "RAD51D archival" "RAD51D relapse"  "BRIP1 archival"  "BRIP1 relapse"
# [9] "PALB2 archival"  "PALB2 relapse"   "BARD1 archival"  "BARD1 relapse"
#[13] "CDK12 archival"  "CDK12 relapse"   "TP53 archival"   "TP53 relapse"
#[17] "KRAS archival"   "KRAS relapse"    "PIK3CA archival" "PIK3CA relapse"
#[21] "NF1 archival"    "NF1 relapse"     "RB1 archival"    "RB1 relapse"
#[25] "NRAS archival"   "NRAS relapse"

	# map short mutation classes to colour labels
	col = c(
	  'shared mutation'= 'cyan',
	  'archival only mutation'= 'purple',
	  'relapse only mutation'= 'gold'
	)
	alter_fun = list(
		background = ComplexHeatmap::alter_graphic("rect", fill = "#CCCCCC"),
		`shared mutation` = ComplexHeatmap::alter_graphic("rect", fill = col["shared mutation"]),
		`archival only mutation` = ComplexHeatmap::alter_graphic("rect", fill = col["archival only mutation"]),
		`relapse only mutation` = ComplexHeatmap::alter_graphic("rect", fill = col["relapse only mutation"])
	)


	print(heatmap_table)
	print(new_somatic_variants %>% colnames()) 
	#quit()

	combined_oncoprint = ComplexHeatmap::oncoPrint(
	  mat=new_somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum resistant'],
	  alter_fun = alter_fun,
	  show_heatmap_legend=FALSE,
	  row_labels=rep('', length(gene_list)),
	  column_title='Pt resistant',
          column_title_gp = gpar(fontsize = column_title_font_size, fontface='bold'),
	  pct_gp = gpar(fontsize = row_label_font_size),
	  row_names_gp = gpar(fontsize=row_label_font_size, fontface='bold'),
	  column_order = colnames(new_somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum resistant']),
	  row_order = gene_list
	) +
	ComplexHeatmap::oncoPrint(
	  mat=new_somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum sensitive'],
	  alter_fun = alter_fun,
	  column_title='Pt sensitive',
          column_title_gp = gpar(fontsize = column_title_font_size, fontface='bold'),
	  pct_gp = gpar(fontsize = row_label_font_size),
	  row_names_gp = gpar(fontsize=row_label_font_size, fontface='bold'),
	  column_order = colnames(new_somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum sensitive']),
	  heatmap_legend_param = list(title='', labels_gp = gpar(fontsize = 6, fontface='bold'), ncol=5, ncol=2, by_row=TRUE),
	  row_order = gene_list
	) 

	combined_oncoprint = ComplexHeatmap::draw(
		combined_oncoprint, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
	)
	
	dev.off()
	return()
}

#generate_somatic_oncoprint('somatic_variants_for_oncorprint.tsv', 'HRD_somatic_oncoprint.pdf', c('TP53','BRCA1','BRCA2','FANCM','BARD1','RAD51B','RAD51C','RAD51D','BRIP1','PALB2'))
generate_somatic_oncoprint(
	somatic_variants=snakemake@input[['data_for_somatic_oncoprint']],
	somatic_oncoprint_output_file=snakemake@output[['somatic_oncoprint']],
	gene_set_analysed=snakemake@params$gene_set_analysed
)
