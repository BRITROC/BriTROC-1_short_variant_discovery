library(magrittr)
library(DBI)
library(RPostgres)
library(ComplexHeatmap)

generate_germline_oncoprint = function(filtered_germline_variants, non_hgsoc_samples, samples_with_no_good_sequencing, samples_with_very_low_purity, analysis_type, output_file) {

	source('~/.Renviron')
	source('functions.R')

	# establish database connections
	britroc_con = make_connection_to_postgres_server('britroc1', 'jblab-db.cri.camres.org', 5432)
	clarity_con = make_connection_to_postgres_server('clarity', 'jblab-db.cri.camres.org', 5432)
	
	# read in germline variant information and reformat information
	germline_variants = readr::read_tsv(filtered_germline_variants) %>% 
		dplyr::rename(patient_id='fk_britroc_number', Consequence='variant_type', SYMBOL='gene_symbol') %>% 
		dplyr::select(patient_id,Consequence,SYMBOL)
	
	# recode germline variant information
	germline_variants$Consequence = dplyr::recode(
		germline_variants$Consequence,
		'splice_acceptor'='splice_acceptor_variant',
		'frameshift'='frameshift_variant',
		'nonsynonymous'='missense_variant',
		'stop_gained,splice_region'='stop_gained,splice_region_variant',
		'splice_region,intron'='splice_region_variant,intron_variant',
		'framehift,splice_region'='frameshift_variant'
	)
	
	# factorise variant consequence column
	germline_variants$Consequence = factor(germline_variants$Consequence,
	                                       levels=c(
	                                         'frameshift_variant',
	                                         'inframe_deletion',
	                                         'inframe_insertion',
	                                         'stop_gained',
	                                         'stop_gained,frameshift_variant',
	                                         'stop_gained,splice_region_variant',
	                                         'missense_variant',
	                                         'splice_donor_variant',
	                                         'splice_acceptor_variant',
	                                         'splice_region_variant,intron_variant',
						 'splice_region_variant,synonymous_variant'
	                                         ))
	
	# reformat column
	germline_variants = germline_variants %>% dplyr::arrange(patient_id,Consequence)
	
	# remove samples that fail QC filters
	relevant_samples = remove_non_relevant_samples(non_hgsoc_samples, samples_with_no_good_sequencing, samples_with_very_low_purity, britroc_con, clarity_con, analysis_type)

	# only retain variants in patients with at least one germline, relapse and archival sample
	germline_variants = germline_variants %>% dplyr::filter(patient_id %in% relevant_samples$fk_britroc_number)
	
	# reformat
	germline_variants =
	  germline_variants %>% dplyr::select(patient_id,SYMBOL,Consequence) %>%
	  unique()
	
	# remove synonymous and intron variants
	# TODO: I am not sure if this line does anything - potentially remove
	germline_variants = 
	  germline_variants %>% dplyr::filter(!Consequence %in% c('intron_variant','synonymous_variant')) %>%
	  unique()
	
	# this information is needed so that the oncoprint shows information even for patients without a mutation in any one gene
	# TODO: We should probably wrap this into a function
	somatic_samples_with_no_mutations = dplyr::anti_join(
	  relevant_samples, germline_variants, by=c('fk_britroc_number'='patient_id')
	)
	
	somatic_samples_with_no_mutations = somatic_samples_with_no_mutations %>% dplyr::rename(patient_id='fk_britroc_number')
	
	# reformat
	somatic_samples_with_no_mutations = somatic_samples_with_no_mutations %>%
	  dplyr::select(patient_id) %>% 
	  unique() %>% 
	  dplyr::mutate(variant_type='dummy', gene_symbol='BRCA1') %>% 
	  dplyr::select(patient_id,variant_type,gene_symbol)
	
	somatic_samples_with_mutations = germline_variants %>% dplyr::select(patient_id,SYMBOL,Consequence)
	somatic_samples_with_no_mutations =  somatic_samples_with_no_mutations %>% dplyr::ungroup()
	colnames(somatic_samples_with_mutations) = c('patient_id','gene_symbol','variant_type')
	
	somatic_variants = rbind(somatic_samples_with_mutations, somatic_samples_with_no_mutations)
	
	# remap variant type categories
	somatic_variants$variant_type = dplyr::recode(
	  somatic_variants$variant_type,
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
	
	# for simplicity for now - only one mutation per patient-gene group
	somatic_variants = somatic_variants %>% tibble::as_tibble() %>% 
	  dplyr::group_by(patient_id,gene_symbol) %>%
	  dplyr::filter(dplyr::row_number()==1) %>% dplyr::ungroup()
	
	# convert the data frame to a matrix
	somatic_variants = somatic_variants %>% 
	  tidyr::pivot_wider(names_from=patient_id, values_from=variant_type) %>% 
	  as.matrix()
	
	row.names(somatic_variants) = somatic_variants[,1]
	somatic_variants = somatic_variants[,-1]
	
	# replace NA values
	somatic_variants = tidyr::replace_na(somatic_variants, replace = '')
	
	print(somatic_variants)
	
	# map short mutation classes to labels
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
	  `stop gained` = ComplexHeatmap::alter_graphic("rect", fill = col["stop_gained"]),
	  `inframe indel` = ComplexHeatmap::alter_graphic("rect", fill = col["inframe indel"]),
	 `splice region SNV` = ComplexHeatmap::alter_graphic("rect", fill = col["splice_region_SNV"]),
	  missense = ComplexHeatmap::alter_graphic("rect", fill = col["missense"])
	)
	
	# generate a platinum sensitivity partition
	heatmap_patients = colnames(somatic_variants)
	heatmap_table = tibble::tibble(fk_britroc_number=as.integer(heatmap_patients))
	print(heatmap_table)
	
	patients = dbReadTable(britroc_con, 'patients')
	print(patients)
	
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
	
	genes_with_variants = c('BRCA1','BRCA2','RAD51B','RAD51C','RAD51D','BRIP1')

	# remove the PALB2 row which does not have any variants
	somatic_variants = somatic_variants[1:6,]

	# set the order variants appear in the oncoprint
	variant_type_order = c('frameshift', 'stop gained', 'splice region SNV', 'missense')

	# order tables by archival TP53 variant type	
	somatic_variants = somatic_variants[, order( factor(somatic_variants[1,], levels=variant_type_order))]
	# ensure heatmap table has the same order as somatic variants
	heatmap_table = heatmap_table[match(colnames(somatic_variants) %>% as.integer, heatmap_table$fk_britroc_number),]       	

	print(somatic_variants)

	#pdf(output_file, width=10, height=5)
	png(output_file, width=1200)

	column_title_font_size = 8
	row_label_font_size = 6	

        colours = structure(c('red','orange','yellow','blue'), names = c("frameshift", "stop gained", "splice region SNV", "missense"))
		
	combined_oncoprint = ComplexHeatmap::oncoPrint(
	  mat=somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum resistant'],
	  alter_fun = alter_fun,
#	  col = colours,
	  show_heatmap_legend=FALSE,
	  row_labels=c('','','','','',''),
	  column_title='Pt resistant',
	  column_title_gp = gpar(fontsize = column_title_font_size, fontface='bold'),
	  pct_gp = gpar(fontsize = row_label_font_size),
	  row_names_gp = gpar(fontsize=row_label_font_size, fontface='bold'),
	  heatmap_legend_param = list(title='', labels_gp = gpar(fontsize = 6, fontface='bold'), ncol=5),
	  #column_order = colnames(somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum resistant']),
	  row_order = genes_with_variants
	) +
	ComplexHeatmap::oncoPrint(
	  mat=somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum sensitive'],
	  alter_fun = alter_fun,
#	  col = colours,
	  column_title='Pt sensitive',
	  show_heatmap_legend=TRUE,
          column_title_gp = gpar(fontsize = column_title_font_size, fontface='bold'),
	  pct_gp = gpar(fontsize = row_label_font_size),
          row_names_gp = gpar(fontsize=row_label_font_size, fontface='bold'),
          heatmap_legend_param = list(title='', labels_gp = gpar(fontsize = 6, fontface='bold'), ncol=5),
          #column_order = colnames(somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum sensitive']),
	  row_order = genes_with_variants
	) 
	#lgd = Legend(labels = 'missense', title = "", legend_gp = gpar(fill = 4))

	combined_oncoprint = ComplexHeatmap::draw(
		combined_oncoprint, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend=TRUE
	)

	#ComplexHeatmap::draw(lgd)
	
	dev.off()	
}

generate_germline_oncoprint(
	filtered_germline_variants=snakemake@input[['filtered_germline_variants']],
        non_hgsoc_samples = snakemake@config[['non_hgsoc_samples']],
        samples_with_no_good_sequencing = snakemake@config[['samples_with_no_good_sequencing']],
        samples_with_very_low_purity = snakemake@config[['samples_with_very_low_purity']],
	analysis_type = 'germline',
	output_file = snakemake@output[['germline_oncoprint']]
)





