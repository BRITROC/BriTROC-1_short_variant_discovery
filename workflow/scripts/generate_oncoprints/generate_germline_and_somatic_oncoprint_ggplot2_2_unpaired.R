library(magrittr)
library(ggplot2)
library(patchwork)

generate_somatic_oncoprint = function(somatic_variants, somatic_oncoprint_output_file, gene_set_analysed) {

	source('~/.Renviron')
	source('functions.R')

	# read in preprocessed data
	somatic_variants = readr::read_tsv(somatic_variants)

	print(somatic_variants)

	# all genes examined in this analysis
	all_genes =  tibble::tibble(names=gene_set_analysed) %>% dplyr::filter(names %in% somatic_variants$gene_symbol)	

	# create a genes column
	somatic_variants_tmp = somatic_variants %>% 
		dplyr::mutate(gene_symbol=stringr::str_extract(gene_symbol, pattern='[A-Z0-9]+_') %>% 
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

	
	# fill in blanks for those patient-gene combinations without a mutation
	# add in data from patients with no mutations in any gene type
	#metadata = readr::read_tsv('config/matched_somatic_metadata.tsv')
	#metadata$fk_britroc_number %>% unique() %>% length() %>% print()

	#non_hgsoc_samples = c('JBLAB-4114','JBLAB-4916','IM_249','IM_250','IM_234','IM_235','IM_236','IM_237','JBLAB-4271','IM_420','IM_262','JBLAB-4922','JBLAB-4923','IM_303','IM_290','IM_43','IM_293','IM_307','IM_308','IM_309','IM_424','IM_302','IM_303','IM_304','IM_305','JBLAB-19320','IM_61','IM_62','IM_63','IM_397','IM_302','IM_98','JBLAB-4210','IM_147','JBLAB-4216','IM_44')
	#samples_with_no_good_sequencing = c('IM_144','IM_435','IM_436','IM_158','IM_296','IM_373','IM_154','IM_297','IM_365','IM_432','IM_429','IM_368','IM_441')
	#samples_with_very_low_purity = c('IM_1','IM_2','IM_3','IM_4','IM_20','IM_26','IM_27','IM_69','IM_86','IM_90','IM_93','IM_94','IM_173','IM_177','IM_179','IM_200','IM_241','IM_242','IM_417','IM_418','IM_419','IM_420','IM_221','IM_264','IM_329','IM_289','IM_308','IM_309','IM_338','IM_339','IM_340','IM_341','IM_342','IM_432','IM_372','IM_272','IM_392','IM_114','IM_382','JBLAB-4911')
	
	#metadata = metadata %>% dplyr::filter(!fk_sample %in% non_hgsoc_samples)
	#metadata = metadata %>% dplyr::filter(!fk_sample %in% samples_with_no_good_sequencing)
	#metadata = metadata %>% dplyr::filter(!fk_sample %in% samples_with_very_low_purity)

	#setdiff(metadata$fk_britroc_number, somatic_variants$patient_id) %>% print()
	#setdiff(somatic_variants$patient_id, metadata$fk_britroc_number) %>% print()
	
	#somatic_variants = somatic_variants %>% dplyr::filter(patient_id %in% metadata$fk_britroc_number)

	#no_mutation_table = tibble::tibble(
	#	patient_id=somatic_variants$patient_id %>% unique(),
	#	variant_type='no mutation'	
	#)

	somatic_variants = somatic_variants %>% tibble::add_row(
		patient_id=191,
		variant_type='no mutation',
		gene_symbol='TP53'
	)

	somatic_variants = somatic_variants %>% tibble::add_row(
		patient_id=191,
		variant_type='no mutation',
		gene_symbol='BARD1'
	)
	
	somatic_variants = somatic_variants %>% dplyr::mutate(variant_class='somatic')
	
	germline_variants = readr::read_tsv(snakemake@input[['germline_data']])
	print(germline_variants)
	germline_variants = germline_variants %>% dplyr::filter(fk_britroc_number %in% somatic_variants$patient_id)
	germline_variants$HGVS_united = dplyr::coalesce(germline_variants$HGVSp, germline_variants$HGVSc)
	germline_variants %>% dplyr::select(SYMBOL,HGVSc,HGVSp) %>% unique() %>% print()
	readr::write_tsv(germline_variants %>% dplyr::select(HGVS_united) %>% unique(), 'tmp.txt')


	germline_variants$HGVS_united = germline_variants$HGVS_united %>% stringr::str_remove('ENS[A-Z0-9]+\\.[0-9]+:')
	germline_variants %>% dplyr::filter(SYMBOL=='BRCA1') %>% .$HGVS_united %>% unique() %>% print()	
	
	# MTBP job code: octopus germline matched 
	germline_variants = 
		germline_variants %>%
		dplyr::filter(
			HGVS_united %in% 
			c(
				'p.Tyr1646Ter',
				'p.Glu23ValfsTer17',
				'p.Val299ArgfsTer4',
				'p.Lys583AsnfsTer3',
				'p.His692MetfsTer9',
				'p.Glu881Ter',
				'p.Lys894ThrfsTer8',
				'p.Asn1002ThrfsTer22',
				'p.Gly1077AlafsTer8',
				'p.Arg1203Ter',
				'p.Ile1486Ter',
				'c.4548-2A>G',
				'p.Trp1529Ter',
				'p.Ala1644Gly',
				'c.302-2del',
				'c.5049+6T>G',
				'p.Val1757Ala',
				'p.Ala1729Glu',
				'p.Glu1035Ter',
				'p.Val220IlefsTer4',
				'p.Pro9GlnfsTer16',
				'p.Gln397LeufsTer25',
				'p.Asp1469LysfsTer11',
				'p.Leu2092ProfsTer7',
				'p.Glu1493ValfsTer10',
				'p.Thr3033LeufsTer29',
				'p.Gln2829Ter',
				'p.Glu1571GlyfsTer3',
				'p.Asn1784HisfsTer2',
				'p.Arg798Ter',
				'p.Asn1039IlefsTer2',
				'p.Arg47Ter',
				'p.Gln133Ter',
				'p.Asp94IlefsTer7'
			)
		)
	germline_variants %>% dplyr::filter(SYMBOL=='BRCA1') %>% .$HGVS_united %>% unique() %>% print()
        germline_variants %>% .$HGVS_united %>% unique() %>% print()
	#quit()
	germline_variants = germline_variants %>% dplyr::select(fk_britroc_number, SYMBOL, Consequence) %>% unique()

	germline_variants$Consequence %>% as.factor() %>% levels() %>% print()
	#quit()

	germline_variants$Consequence = germline_variants$Consequence %>%
		dplyr::recode(
			'frameshift_variant'='frameshift',
			'frameshift_variant,splice_region_variant'='frameshift',
			'inframe_deletion'='inframe deletion',
			'missense_variant'='missense',
			'missense_variant,splice_region_variant'='missense',
			'splice_region_variant,intron_variant'='splice region SNV',
			'stop_gained,splice_region_variant'='stop gained',
			'splice_acceptor'='splice region SNV',
			'splice_acceptor_variant'='splice region SNV',
			'stop_gained'='stop gained',
			'stop_gained,splice_region'='stop gained'
		)

	germline_variants$Consequence %>% as.factor() %>% levels() %>% print()
	#quit()

	germline_variants = germline_variants %>% dplyr::mutate(variant_class='germline')
	germline_variants = germline_variants %>% dplyr::rename(patient_id=fk_britroc_number, gene_symbol=SYMBOL, variant_type=Consequence)

	print(somatic_variants)
	print(germline_variants)

	somatic_variants = rbind(somatic_variants,germline_variants)
	print(somatic_variants)
	#quit()

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

	no_mutation_table = no_mutation_table %>% dplyr::mutate(variant_class='somatic')

	# combine both variant type tables together
	somatic_variants = rbind(
		somatic_variants,
		no_mutation_table
	) %>% unique()

	print(germline_variants, n=Inf)
	#quit()

	# enforce variant type legend order and gene order
	somatic_variants$variant_type = factor(
		somatic_variants$variant_type,
		levels=c('frameshift','stop gained','splice region SNV','missense','inframe indel','no mutation')
	)
	
	somatic_variants$gene_symbol = factor(
		somatic_variants$gene_symbol,
		levels=rev(
			c('TP53','BRCA1','BRCA2','BRIP1','BARD1','PALB2','RAD51B','RAD51C'))
	)

	#quit()

	# convert patient_id to a factor
	somatic_variants$patient_id = somatic_variants$patient_id %>% as.factor()
	somatic_variants$gene_symbol = somatic_variants$gene_symbol %>% as.factor()

	somatic_variants$variant_type_TP53 = 
		dplyr::if_else(
			somatic_variants$gene_symbol %in% c('TP53'),
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
		somatic_variants$gene_symbol %in% c('BRCA1'),
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
			somatic_variants$gene_symbol %in% c('BRCA2'),
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

	somatic_variants %>% print(width=Inf)
	#quit()

	gene_percentages = 
		somatic_variants %>%
		dplyr::filter(variant_type!='no mutation') %>%
		dplyr::group_by(gene_symbol) %>%
		dplyr::select(gene_symbol) %>% 
		table(gene_symbol = .) %>%
		data.frame() %>%
		#dplyr::summarise(n = dplyr::n()) %>%
		dplyr::mutate(percentage = Freq / length(somatic_variants$patient_id %>% unique())) %>%
		.$percentage %>%
                `*`(100) %>%  
		round(digits=0) %>%
		as.character() %>%
		paste('%',sep='')

	print(gene_percentages)
	class(gene_percentages) %>% print()	

	somatic_variants %>% dplyr::filter(grepl('BRCA1',gene_symbol)) %>% print()
	#quit()

	somatic_variants = somatic_variants %>% dplyr::group_by(patient_id, gene_symbol) %>%
		dplyr::slice_head(n=1) %>% dplyr::ungroup()

	somatic_variants$alpha_value = dplyr::if_else(
		somatic_variants$variant_class=='somatic',
		1.0,
		0.5
	)

	p1 = ggplot(somatic_variants, aes(x=patient_id, y=as.numeric(gene_symbol))) +
                geom_tile(aes(fill = variant_type), colour='white', size=0.9, alpha=somatic_variants$alpha_value) +
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
			breaks = 1:length(somatic_variants$gene_symbol %>% as.factor() %>% levels()),
			labels = 
				somatic_variants$gene_symbol %>% as.factor() %>% levels(),
				sec.axis = sec_axis(
					~.,
					breaks = 1:length(somatic_variants$gene_symbol %>% as.factor %>% levels()),
					labels = gene_percentages
					)
		)

	patient_table = readr:read_tsv(snakemake@input['patient_table_file_path'])
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

	p_final = p_final + plot_annotation(
		title = 'BriTROC-1: octopus somatic (bold) + octopus germline (pale) variants',
		subtitle = 'Calling mode: Matched (non-tumour blood sample) and unpaired analysis'
	)

	ggsave(somatic_oncoprint_output_file, p_final, dev='png', width=16.0, height=9.0, scale=0.75, dpi=300)

	return()
}

generate_somatic_oncoprint(
	somatic_variants=snakemake@input[['data_for_somatic_oncoprint']],
	somatic_oncoprint_output_file=snakemake@output[[1]],
	gene_set_analysed=snakemake@params$gene_set_analysed
)
