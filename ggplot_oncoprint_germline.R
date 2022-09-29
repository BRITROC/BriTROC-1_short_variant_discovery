# A ggplot implementation of oncoprints
# this time for germline variants

library(patchwork)
library(ggplot2)

#theme_set(theme_gray(base_size = 19))

# read in data

oncoprint_data = readr::read_tsv('results/variant_analysis/non_TP53/panel_6_28/collated/final_germline_variants.tsv')

print(oncoprint_data)
oncoprint_data %>% colnames() %>% print()

# remove NA results
oncoprint_data = oncoprint_data %>% dplyr::filter(!is.na(variant_type))

oncoprint_data$gene_symbol %>% as.factor %>% levels()

gene_list = c("BRCA1","BRCA2","BRIP1","PALB2","RAD51B","RAD51C","RAD51D")

print(oncoprint_data)

# order table
oncoprint_data = oncoprint_data %>% dplyr::arrange(fk_britroc_number, gene_symbol)

# convert britroc number to a factor
oncoprint_data$britroc_number = as.factor(oncoprint_data$fk_britroc_number)

# transform data into a single value per patient-gene combination
#transform_data = function(patient, gene) {
#	combined_table_tmp = combined_table %>% dplyr::filter(grepl(gene, gene_list)) %>% dplyr::filter(britroc_number==patient)
#
#	new_combined_table_tmp = tibble::tibble(britroc_number=patient, gene_list=gene, variant_type=NA)
#
#	new_combined_table_tmp$variant_type  = dplyr::if_else(combined_table_tmp$variant_type[1] == 'mutation_present' & combined_table_tmp$variant_type[2] == 'mutation_present', 'shared mutation' , as.character(NA))
#	new_combined_table_tmp$variant_type  = dplyr::if_else(combined_table_tmp$variant_type[1] == 'mutation_present' & combined_table_tmp$variant_type[2] == 'mutation_not_present', 'archival only mutation' , new_combined_table_tmp$variant_type)
#	new_combined_table_tmp$variant_type  = dplyr::if_else(combined_table_tmp$variant_type[1] == 'mutation_not_present' & combined_table_tmp$variant_type[2] == 'mutation_present', 'relapse only mutation' , new_combined_table_tmp$variant_type)
#	new_combined_table_tmp$variant_type  = dplyr::if_else(combined_table_tmp$variant_type[1] == 'mutation_not_present' & combined_table_tmp$variant_type[2] == 'mutation_not_present', 'no mutation' , new_combined_table_tmp$variant_type)
#
#	return(new_combined_table_tmp)
#
#}

#combined_table_reduced = combined_table %>% dplyr::mutate(gene_list = gene_list %>% stringr::str_remove('_archival') %>% stringr::str_remove('_relapse')) %>%
#	dplyr::select(-variant_type) %>% unique()

# remove all genes which only have no_mutations
#mutation_freq = new_combined_table %>% dplyr::group_by(gene_list,variant_type) %>% dplyr::summarise(n=dplyr::n()) %>% dplyr::group_by(gene_list) %>%
#	dplyr::summarise(n=dplyr::n()) %>% dplyr::filter(n==1)
#print(mutation_freq, n=Inf)
#new_combined_table = new_combined_table %>% dplyr::filter(!gene_list %in% mutation_freq$gene_list)

# factories gene list
#new_combined_table$gene_list = factor(new_combined_table$gene_list, 
#	levels=c('TP53','BRCA1','BRCA2','BRIP1','RAD51D','PALB2','BARD1','CDK12','KRAS','PIK3CA','NF1','RB1','NRAS') %>% rev()
#	)


#new_combined_table2 = dplyr::inner_join(new_combined_table, patient_table, by='britroc_number')
#new_combined_table2$variant_type = factor(new_combined_table2$variant_type, levels=c('shared mutation','archival only mutation','relapse only mutation','no mutation'))

# strictly for the purpose of matrix ordering and nothing else
#new_combined_table2$variant_type_no_tp53 = new_combined_table2$variant_type

#print(new_combined_table2) 

#new_combined_table2$gene_list_not_factor = new_combined_table2$gene_list %>% as.character()

#print(new_combined_table2) 

# create table

# ensure the britroc number is not being used as a factor 
#new_combined_table2$britroc_number = new_combined_table2$britroc_number %>% as.integer()

#new_combined_table2 = new_combined_table2 %>% dplyr::arrange(variant_type_BRCA1)

gene_list=oncoprint_data$gene_symbol %>% as.factor %>% levels()

oncoprint_data = oncoprint_data %>% dplyr::select(fk_britroc_number,gene_symbol,variant_type)

print(oncoprint_data)

# make connection to database
library(DBI)
library(RPostgres)

readRenviron('~/.Renviron')
source('functions.R')

britroc_con <- dbConnect(RPostgres::Postgres(),
                 dbname='britroc1',
                 host='jblab-db.cri.camres.org',
                 port = 5432,
                 user = Sys.getenv('jblab_db_username'),
                 password = Sys.getenv('jblab_db_password')
)

clarity_con <- dbConnect(RPostgres::Postgres(),
                 dbname='britroc1',
                 host='jblab-db.cri.camres.org',
                 port = 5432,
                 user = Sys.getenv('jblab_db_username'),
                 password = Sys.getenv('jblab_db_password')
)

non_hgsoc_samples = c('JBLAB-4114','JBLAB-4916','IM_249','IM_250','IM_234','IM_235','IM_236','IM_237','JBLAB-4271','IM_420','IM_262','JBLAB-4922','JBLAB-4923','IM_303','IM_290','IM_43','IM_293','IM_307','IM_308','IM_309','IM_424','IM_302','IM_303','IM_304','IM_305','JBLAB-19320','IM_61','IM_62','IM_63','IM_397','IM_302','IM_98','JBLAB-4210','IM_147','JBLAB-4216','IM_44')
samples_with_no_good_sequencing = c('IM_144','IM_435','IM_436','IM_158','IM_296','IM_373','IM_154','IM_297','IM_365','IM_432','IM_429','IM_368','IM_441')
samples_with_very_low_purity = c('IM_1','IM_2','IM_3','IM_4','IM_20','IM_26','IM_27','IM_69','IM_86','IM_90','IM_93','IM_94','IM_173','IM_177','IM_179','IM_200','IM_241','IM_242','IM_417','IM_418','IM_419','IM_420','IM_221','IM_264','IM_329','IM_289','IM_308','IM_309','IM_338','IM_339','IM_340','IM_341','IM_342','IM_432','IM_372','IM_272','IM_392','IM_114','IM_382','JBLAB-4911')

# remove samples that fail QC filters
relevant_samples = remove_non_relevant_samples(non_hgsoc_samples, samples_with_no_good_sequencing, samples_with_very_low_purity, britroc_con, clarity_con, 'germline')

# only include patients with relevant samples
samples = dbReadTable(britroc_con, 'sample') %>% tibble::as_tibble()
samples = samples %>% dplyr::filter(type=='germline')

dim(samples)
samples = samples %>% dplyr::filter(fk_britroc_number %in% relevant_samples$fk_britroc_number)
dim(samples)

patients_with_germline_samples = samples$fk_britroc_number %>% as.factor %>% levels()

print(patients_with_germline_samples)

print(oncoprint_data)

patients_without_a_germline_mutation = tibble::tibble(
	fk_britroc_number=patients_with_germline_samples, variant_type='no mutation')
gene_table = tibble::tibble(
	gene_symbol=oncoprint_data$gene_symbol %>% as.factor() %>% levels()
	)

patients_without_a_germline_mutation = 
	dplyr::inner_join(
		patients_without_a_germline_mutation,
		gene_table,
		by = character()
	)

oncoprint_data$fk_britroc_number = oncoprint_data$fk_britroc_number %>% as.factor()
patients_without_a_germline_mutation$fk_britroc_number = 
	patients_without_a_germline_mutation$fk_britroc_number %>% as.factor()

patients_without_a_germline_mutation = 
	dplyr::anti_join(
		patients_without_a_germline_mutation,
		oncoprint_data,
		by=c('fk_britroc_number','gene_symbol')
	)

oncoprint_data = rbind(
		oncoprint_data,
		patients_without_a_germline_mutation
	)

oncoprint_data$variant_type = oncoprint_data$variant_type %>% 
	dplyr::recode(
		'frameshift,splice_region'='frameshift',
		'nonsynonymous'='missense',
		'splice_acceptor'='splice region SNV',
		'splice_region,intron'='splice region SNV',
		'stop_gained,splice_region'='splice region SNV',
		'stop_gained'='stop gained'	
	)

oncoprint_data$variant_type %>% as.factor %>% levels %>% print()

oncoprint_data$variant_type = factor(oncoprint_data$variant_type, levels=c('frameshift','stop gained','splice region SNV','missense','no mutation'))
oncoprint_data$gene_symbol = factor(oncoprint_data$gene_symbol, levels=rev(c('BRCA1','BRCA2','BRIP1','PALB2','RAD51B','RAD51C','RAD51D')))

oncoprint_data$variant_type_BRCA1 = 
	dplyr::if_else(
	oncoprint_data$gene_symbol == 'BRCA1',
	oncoprint_data$variant_type %>% as.character(),
	'no mutation'
)

oncoprint_data$variant_type_BRCA1 = oncoprint_data$variant_type_BRCA1 %>%
	dplyr::recode(
	'frameshift'=1L,
	'missense'=4L,
	'splice region SNV'=3L,
	'stop gained'=2L,
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
	'frameshift'=5L,
	'missense'=8L,
	'splice region SNV'=7L,
	'stop gained'=6L,
	'no mutation'=200L
	)

oncoprint_data$variant_type_BRIP1 = 
	dplyr::if_else(
	oncoprint_data$gene_symbol == 'BRIP1',
	oncoprint_data$variant_type %>% as.character(),
	'no mutation'
)

oncoprint_data$variant_type_BRIP1 = oncoprint_data$variant_type_BRIP1 %>%
	dplyr::recode(
	'frameshift'=9L,
	'missense'=12L,
	'splice region SNV'=11L,
	'stop gained'=10L,
	'no mutation'=200L
	)

oncoprint_data$variant_type_BRCA1_BRCA2_BRIP1 = 
	oncoprint_data$variant_type_BRCA1 + 
	oncoprint_data$variant_type_BRCA2 +
	oncoprint_data$variant_type_BRIP1

oncoprint_data$fk_britroc_number = 
	forcats::fct_reorder(
		oncoprint_data$fk_britroc_number, 
		oncoprint_data$variant_type_BRCA1_BRCA2_BRIP1, 
		.fun=min
)

patient_table = dbReadTable(britroc_con, 'patients') %>% tibble::as_tibble()
patient_table$britroc_number = as.factor(patient_table$britroc_number)

oncoprint_data = dplyr::inner_join(oncoprint_data, patient_table, by=c('fk_britroc_number'='britroc_number'))
oncoprint_data$pt_sensitivity_at_reg = oncoprint_data$pt_sensitivity_at_reg %>% as.factor()

new_combined_table = oncoprint_data %>% dplyr::select(fk_britroc_number,pt_sensitivity_at_reg,variant_type_BRCA1_BRCA2_BRIP1) %>% 
	dplyr::mutate(goo='Histological type') %>% unique()

new_combined_table$fk_britroc_number = 
	forcats::fct_reorder(
		new_combined_table$fk_britroc_number, 
		new_combined_table$variant_type_BRCA1_BRCA2_BRIP1, 
		.fun=min
)


p1 = ggplot(oncoprint_data, aes(x=fk_britroc_number, y=gene_symbol)) + 
		geom_tile(aes(fill = variant_type), colour = "white", size=0.9) +
		scale_fill_manual(values=c('#604187','#FF1493','#FFD700','#00FFFF','#D3D3D3')) +
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
		) + ylab('')

p5 = ggplot(new_combined_table, 
		aes(x=fk_britroc_number, y=goo, width=1.0)) +  geom_tile(aes(fill = pt_sensitivity_at_reg), colour = "white") +
		scale_fill_manual(values=c('#DCE319FF', '#20A387FF')) +
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

ggsave('BriTROC-1_manuscript2_germline_oncoprint.png', p_final, device='png', width=6.5, height=4.0, scale=1.50, dpi=300) # height=5.4
