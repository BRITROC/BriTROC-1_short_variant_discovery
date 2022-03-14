# A ggplot implementation of oncoprints

library(patchwork)
library(ggplot2)

theme_set(theme_gray(base_size = 22))

# read in data
oncoprint_data = readr::read_tsv('results/data_for_whole_cohort_oncoprint_intercalated.targeted.tsv')

# remove NA results
oncoprint_data = oncoprint_data %>% dplyr::filter(!is.na(variant_type))

metadata_sheet = readr::read_tsv('config/tumour_metadata_with_one_of_both_types.tsv')

gene_list = c("BRCA1","BRCA2","RAD51C","RAD51D","RAD51B","BRIP1","FANCM","PALB2","BARD1","CDK12","EGFR","PTEN","TP53","KRAS","BRAF","PIK3CA","CTNNB1","NF1","RB1","NRAS")

gene_list = c(
	paste(gene_list, 'archival', sep='_'),
	paste(gene_list, 'relapse', sep='_')
)

oncoprint_data$variant_type = 'mutation_present'

# make connection to database
library(DBI)
library(RPostgres)

readRenviron('~/.Renviron')
britroc_con <- dbConnect(RPostgres::Postgres(),
                 dbname='britroc1',
                 host='jblab-db.cri.camres.org',
                 port = 5432,
                 user = Sys.getenv('jblab_db_username'),
                 password = Sys.getenv('jblab_db_password')
)

patient_table = dbReadTable(britroc_con, 'patients') %>% tibble::as_tibble()
patient_table = patient_table %>% dplyr::select(britroc_number, age, histological_type, tumour_location_at_diagnosis, tumour_stage_at_diagnosis, 
					pt_sensitivity_at_reg)

patient_table$britroc_number = as.factor(patient_table$britroc_number)

#patient_table$tumour_location_at_diagnosis %>% table()
#quit()

patient_table$tumour_stage_at_diagnosis = patient_table$tumour_stage_at_diagnosis %>% as.factor()

print(oncoprint_data)

print(metadata_sheet$fk_britroc_number %>% as.factor %>% levels)

print(gene_list)

# create a cross-join between two vectors
gene_list_tibble = tibble::tibble(gene_list=gene_list)
britroc_number_tibble = tibble::tibble(britroc_number = metadata_sheet$fk_britroc_number %>% as.factor %>% levels %>% as.numeric())

mutation_not_present_table = dplyr::full_join(gene_list_tibble, britroc_number_tibble, by=character())
mutation_not_present_table = mutation_not_present_table %>% dplyr::mutate(variant_type='mutation_not_present')

print(mutation_not_present_table)

mutation_not_present_table = dplyr::anti_join(mutation_not_present_table, oncoprint_data, by=c('britroc_number'='patient_id','gene_list'='gene_symbol_tumour_type'))

print(mutation_not_present_table)

# combine the mutation_present and the mutation_non_present data tables
oncoprint_data = oncoprint_data %>% dplyr::rename(gene_list='gene_symbol_tumour_type', britroc_number='patient_id')

print(oncoprint_data)

combined_table = rbind(oncoprint_data, mutation_not_present_table)

# order table
combined_table = combined_table %>% dplyr::arrange(britroc_number, gene_list)

# convert britroc number to a factor
combined_table$britroc_number = as.factor(combined_table$britroc_number)

# transform data into a single value per patient-gene combination
transform_data = function(patient, gene) {
	combined_table_tmp = combined_table %>% dplyr::filter(grepl(gene, gene_list)) %>% dplyr::filter(britroc_number==patient)

	new_combined_table_tmp = tibble::tibble(britroc_number=patient, gene_list=gene, variant_type=NA)

	new_combined_table_tmp$variant_type  = dplyr::if_else(combined_table_tmp$variant_type[1] == 'mutation_present' & combined_table_tmp$variant_type[2] == 'mutation_present', 'shared mutation' , as.character(NA))
	new_combined_table_tmp$variant_type  = dplyr::if_else(combined_table_tmp$variant_type[1] == 'mutation_present' & combined_table_tmp$variant_type[2] == 'mutation_not_present', 'archival only mutation' , new_combined_table_tmp$variant_type)
	new_combined_table_tmp$variant_type  = dplyr::if_else(combined_table_tmp$variant_type[1] == 'mutation_not_present' & combined_table_tmp$variant_type[2] == 'mutation_present', 'relapse only mutation' , new_combined_table_tmp$variant_type)
	new_combined_table_tmp$variant_type  = dplyr::if_else(combined_table_tmp$variant_type[1] == 'mutation_not_present' & combined_table_tmp$variant_type[2] == 'mutation_not_present', 'no mutation' , new_combined_table_tmp$variant_type)

	return(new_combined_table_tmp)

}

combined_table_reduced = combined_table %>% dplyr::mutate(gene_list = gene_list %>% stringr::str_remove('_archival') %>% stringr::str_remove('_relapse')) %>%
	dplyr::select(-variant_type) %>% unique()

new_combined_table=purrr::map2_dfr(.x=combined_table_reduced$britroc_number, .y=combined_table_reduced$gene_list, .f=transform_data)

print(new_combined_table)

# remove all genes which only have no_mutations
mutation_freq = new_combined_table %>% dplyr::group_by(gene_list,variant_type) %>% dplyr::summarise(n=dplyr::n()) %>% dplyr::group_by(gene_list) %>%
	dplyr::summarise(n=dplyr::n()) %>% dplyr::filter(n==1)
print(mutation_freq, n=Inf)
new_combined_table = new_combined_table %>% dplyr::filter(!gene_list %in% mutation_freq$gene_list)

# factories gene list
new_combined_table$gene_list = factor(new_combined_table$gene_list, 
	levels=c('TP53','BRCA1','BRCA2','BRIP1','RAD51D','PALB2','BARD1','CDK12','KRAS','PIK3CA','NF1','RB1','NRAS') %>% rev()
	)


new_combined_table2 = dplyr::inner_join(new_combined_table, patient_table, by='britroc_number')
new_combined_table2$variant_type = factor(new_combined_table2$variant_type, levels=c('shared mutation','archival only mutation','relapse only mutation','no mutation'))
new_combined_table2 = new_combined_table2 %>% dplyr::arrange(variant_type)
#new_combined_table2$britroc_number = factor(new_combined_table2$britroc_number, levels=new_combined_table2$britroc_number %>% unique())

print(new_combined_table2)

quit()

p1 = ggplot(new_combined_table2, aes(x=britroc_number, y=gene_list, width=1.0)) +  geom_tile(aes(fill = variant_type), colour = "white", size=1.3) +
	     scale_fill_manual(values=c('#FFD700','#D3D3D3','#604187','#00FFFF')) +
		theme(
			axis.title.x=element_blank(),
			axis.text.x=element_blank(),
			axis.ticks.x=element_blank(),
			legend.title=element_blank(),
			axis.text.y = element_text(face="bold"), 
			plot.margin=margin(0.1, 0.1, 1.25, 0.1, 'cm'),
			legend.margin=margin(0, 0, 0, 0),
			legend.box.margin=margin(-10,-10,-10,-10)
		) + ylab('')

# integerate additional patient information

new_combined_table = new_combined_table %>% dplyr::select(britroc_number) %>% dplyr::mutate(goo='Histological type') %>% unique()

new_combined_table = dplyr::inner_join(new_combined_table, patient_table, by='britroc_number')

new_combined_table$pt_sensitivity_at_reg = new_combined_table$pt_sensitivity_at_reg %>% as.factor()
new_combined_table$histological_type = new_combined_table$histological_type %>% as.factor()

new_combined_table$histological_type = new_combined_table$histological_type %>%
	dplyr::recode(
	'endometrioid (grade 3)'='HGEC',
	'high grade serous'='HGSC'
	)

new_combined_table$pt_sensitivity_at_reg = new_combined_table$pt_sensitivity_at_reg %>%
	dplyr::recode(
	'resistant'='Pt resistant',
	'sensitive'='Pt sensitive'
	)

new_combined_table = new_combined_table %>% dplyr::arrange(pt_sensitivity_at_reg)
new_combined_table$britroc_number = factor(new_combined_table$britroc_number, levels=new_combined_table$britroc_number)

print(new_combined_table, width=Inf)

p2 = ggplot(new_combined_table, 
		aes(x=britroc_number, y=goo, width=1.0)) +  geom_tile(aes(fill = histological_type), colour = "white") +
		scale_fill_manual(values=c('#088da5','#d8ce8d')) +
		theme(
			axis.title.x=element_blank(),
			axis.text.x=element_blank(),
			axis.ticks.x=element_blank(),
			legend.title=element_blank(),
			legend.direction='horizontal',
			axis.text.y = element_text(face="bold"),
			plot.margin=margin(0.001, 0.001, 0.001, 0.001, 'cm'),
			legend.margin=margin(0, 0, 0, 0),
			legend.box.margin=margin(-10,-10,-10,-10)
		) + ylab('')

new_combined_table$goo = new_combined_table$goo %>% dplyr::recode('Histological type'='Age')

p3 = ggplot(new_combined_table, 
		aes(x=britroc_number, y=goo, width=1.0)) +  geom_tile(aes(fill = age), colour = "white") +
		scale_fill_continuous(low='#e9d8f2', high='#6a0dad', labels = c(50,80), breaks=c(50,80)) +
		theme(
			axis.title.x=element_blank(),
			axis.text.x=element_blank(),
			axis.ticks.x=element_blank(),
			legend.direction='horizontal',
			legend.title=element_blank(),
			axis.text.y = element_text(face="bold"),
			plot.margin=margin(0.001, 0.001, 0.001, 0.001, 'cm'),
			legend.margin=margin(0, 0, 0, 0),
			legend.box.margin=margin(-10,-10,-10,-10)
		) + ylab('')

new_combined_table$goo = new_combined_table$goo %>% dplyr::recode('Age'='Tumour stage')


p4 = ggplot(new_combined_table, 
		aes(x=britroc_number, y=goo, width=1.0)) +  geom_tile(aes(fill = tumour_stage_at_diagnosis), colour = "white") +
		scale_fill_manual(values=c('#FFC300','#FF5733','#C70039','#900C3F')) +
		theme(
			axis.title.x=element_blank(),
			axis.text.x=element_blank(),
			axis.ticks.x=element_blank(),
			legend.direction='horizontal',
			legend.title=element_blank(),
			axis.text.y = element_text(face="bold"),
			plot.margin=margin(0.001, 0.001, 0.001, 0.001, 'cm'),
			legend.margin=margin(0, 0, 0, 0),
			legend.box.margin=margin(-10,-10,-10,-10)
		) + ylab('')


p5 = ggplot(new_combined_table, 
		aes(x=britroc_number, y=goo, width=1.0)) +  geom_tile(aes(fill = pt_sensitivity_at_reg), colour = "white") +
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
			legend.box.margin=margin(-10,-10,-10,-10)
		) + ylab('')

chemotherapy_lines = dbReadTable(britroc_con, 'chemotherapy_lines') %>% tibble::as_tibble()
chemotherapy_lines = chemotherapy_lines %>% dplyr::group_by(fk_britroc_number) %>% dplyr::summarise(num_lines=max(chemotherapy_line))

chemotherapy_lines$fk_britroc_number = as.factor(chemotherapy_lines$fk_britroc_number)

print(chemotherapy_lines)

chemotherapy_lines$num_lines = dplyr::if_else(chemotherapy_lines$num_lines>4, 5, as.double(chemotherapy_lines$num_lines))

print(chemotherapy_lines, n=Inf)

chemotherapy_lines$num_lines = as.factor(chemotherapy_lines$num_lines)
chemotherapy_lines$num_lines = chemotherapy_lines$num_lines %>% dplyr::recode('5'='4+')

new_combined_table = dplyr::inner_join(new_combined_table, chemotherapy_lines, by=c('britroc_number'='fk_britroc_number')) 

print(new_combined_table, n=Inf)

new_combined_table$goo = new_combined_table$goo %>% dplyr::recode('Tumour stage'='Lines of chemotherapy')

p6 = ggplot(new_combined_table, 
		aes(x=britroc_number, y=goo, width=1.0)) +  geom_tile(aes(fill = num_lines), colour = "white") +
		scale_fill_manual(values=c('#FFC0CB','#FFB6C1','#FF69B4','#FF1493','#C71585')) +
		theme(
			axis.title.x=element_blank(),
			axis.text.x=element_blank(),
			axis.ticks.x=element_blank(),
			legend.title=element_blank(),
			legend.direction='horizontal',
			axis.text.y = element_text(face="bold"),
			plot.margin=margin(0.001, 0.001, 0.001, 0.001, 'cm'),
			legend.margin=margin(0, 0, 0, 0),
			legend.box.margin=margin(-10,-10,-10,-10)
		) + ylab('')

p_final = p5 / p1 / p2 / p4 / p3 / p6  + plot_layout(heights = c(1, 16, 1, 1, 1, 1))

ggsave('tmp3.png', p_final, device='png', width=33, height=16, scale=0.65)
