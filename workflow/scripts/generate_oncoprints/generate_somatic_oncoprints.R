library(magrittr)
library(DBI)
library(RPostgres)
library(ComplexHeatmap)

generate_somatic_oncoprint = function(somatic_variants, somatic_oncoprint_output_file) {

	source('~/.Renviron')
	source('functions.R')

	# establish database connections
	britroc_con = make_connection_to_postgres_server('britroc1', 'jblab-db.cri.camres.org', 5432)

	# read in preprocessed data
	somatic_variants = readr::read_tsv(somatic_variants)

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
	  'stop_gained'= 'orange',
	  'splice_region_SNV'= 'yellow',
	  'inframe indel'= 'cyan',
	  'missense'= 'blue'
	)
	alter_fun = list(
		background = ComplexHeatmap::alter_graphic("rect", fill = "#CCCCCC"),
		frameshift = ComplexHeatmap::alter_graphic("rect", fill = col["frameshift"]),
		stop_gained = ComplexHeatmap::alter_graphic("rect", fill = col["stop_gained"]),
		`inframe indel` = ComplexHeatmap::alter_graphic("rect", fill = col["inframe indel"]),
		splice_region_SNV = ComplexHeatmap::alter_graphic("rect", fill = col["splice_region_SNV"]),
		missense = ComplexHeatmap::alter_graphic("rect", fill = col["missense"])
	)

	# create a dummy matrix including all genes in which a variant has not been found	
	dummy = matrix('', nrow=11, ncol=112)
	row.names(dummy) = c(
		'BARD1_relapse',
		'RAD51B_archival','RAD51B_relapse',
		'RAD51C_archival','RAD51C_relapse',
		'RAD51D_archival','RAD51D_relapse',
		'BRIP1_archival','BRIP1_relapse',
		'PALB2_archival','PALB2_relapse'
	)

	somatic_variants = rbind(somatic_variants, dummy)
	
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

	somatic_variants <- somatic_variants[, order(as.integer(colnames(somatic_variants)))]
	
	# include all relevant genes on the basis of tumour type
	genes_with_variants = c(
	'TP53_archival','TP53_relapse',
	'BRCA1_archival','BRCA1_relapse',
	'BRCA2_archival','BRCA2_relapse',
	'FANCM_archival','FANCM_relapse',
	'BARD1_archival','BARD1_relapse',
	'RAD51B_archival','RAD51B_relapse',
	'RAD51C_archival','RAD51C_relapse',
	'RAD51D_archival','RAD51D_relapse',
	'BRIP1_archival','BRIP1_relapse',
	'PALB2_archival','PALB2_relapse'
	)

	heatmap_table = heatmap_table %>% dplyr::arrange(fk_britroc_number)

	print(somatic_oncoprint_output_file)
	pdf(somatic_oncoprint_output_file)	
	
	combined_oncoprint = ComplexHeatmap::oncoPrint(
	  mat=somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum resistant'],
	  alter_fun = alter_fun,
	  show_heatmap_legend=FALSE,
	  row_labels=c('','','','','','','','','','','','','','','','','','','',''),
	  column_title='platinum resistant',
	  row_order = genes_with_variants
	) +
	ComplexHeatmap::oncoPrint(
	  mat=somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum sensitive'],
	  alter_fun = alter_fun,
	  column_title='platinum sensitive',
	  #heatmap_legend_param = list(labels=c('frameshift','stop gained','splice region SNV','missense'),
	  #                             title='Alterations'),
	  row_order = genes_with_variants
	) 

	# essential step when the pdf function is used inside of a function 
	print(combined_oncoprint)
	
	dev.off()
	return()
}

#generate_somatic_oncoprint('somatic_variants_for_oncorprint.tsv', 'HRD_somatic_oncoprint.pdf')
generate_somatic_oncoprint(
	somatic_variants=snakemake@input[['data_for_somatic_oncoprint']],
	somatic_oncoprint_output_file=snakemake@output[['somatic_oncoprint']]
)
