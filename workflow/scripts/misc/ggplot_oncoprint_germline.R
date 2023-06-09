# A ggplot implementation of oncoprints
# this time for germline variants

library(patchwork)
library(ggplot2)

# read in data
oncoprint_data = readr::read_tsv('results/variant_analysis/non_TP53/panel_6_28/collated/final_germline_variants.tsv')

# remove NA results
oncoprint_data = oncoprint_data %>% dplyr::filter(!is.na(variant_type))
gene_list = c("BRCA1","BRCA2","BRIP1","PALB2","RAD51B","RAD51C","RAD51D")

# order table
oncoprint_data = oncoprint_data %>% dplyr::arrange(fk_britroc_number, gene_symbol)

# convert britroc number to a factor
oncoprint_data$britroc_number = as.factor(oncoprint_data$fk_britroc_number)
gene_list=oncoprint_data$gene_symbol %>% as.factor %>% levels()

oncoprint_data = oncoprint_data %>% dplyr::select(fk_britroc_number,gene_symbol,variant_type)

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
samples = samples %>% dplyr::filter(fk_britroc_number %in% relevant_samples$fk_britroc_number)

patients_with_germline_samples = samples$fk_britroc_number %>% as.factor %>% levels()

print(patients_with_germline_samples %>% length())

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

print(oncoprint_data)
oncoprint_data$fk_britroc_number %>% as.factor %>% levels() %>% length()

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

gene_symbols = oncoprint_data$gene_symbol %>% levels()

oncoprint_data %>% print()

num_patients = oncoprint_data$fk_britroc_number %>% unique() %>% length()
print(num_patients)

oncoprint_tmp = oncoprint_data %>%
	dplyr::filter(variant_type!='no mutation') %>%
	dplyr::group_by(fk_britroc_number,gene_symbol) %>%
	dplyr::slice_head(n=1)

oncoprint_tmp$variant_type %>% as.factor %>% levels()

oncoprint_tmp %>% .$gene_symbol %>% table()

gene_symbol_percentages=c('11%','12%','1%','1%','1%','1%','1%') %>% rev()

p1 = ggplot(oncoprint_data, aes(x=fk_britroc_number, y=as.numeric(gene_symbol))) + 
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
		) + ylab('') +
		scale_y_continuous(breaks = 1:length(gene_symbols),
                     labels = gene_symbols,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(gene_symbols),
                                         labels = gene_symbol_percentages))

new_combined_table$pt_sensitivity_at_reg = new_combined_table$pt_sensitivity_at_reg %>% 
	dplyr::recode(
		'sensitive'='Pt sensitive',
		'resistant'='Pt resistant'
	)

p5 = ggplot(new_combined_table, 
		aes(x=fk_britroc_number, y=goo, width=1.0)) +  geom_tile(aes(fill = pt_sensitivity_at_reg), colour = "white") +
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

ggsave('BriTROC-1_manuscript2_germline_oncoprint.png', p_final, device='png', width=16.0, height=9.0, scale=0.75, dpi=300) # height=5.4
