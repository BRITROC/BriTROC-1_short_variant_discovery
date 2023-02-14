# A script to collate all amplicon panel information together into one table

amplicon_panel_paths = 
	tibble::tribble(
	~amplicon_panel_id, ~amplicon_panel_path,
	'1',	'/scratcha/jblab/amplicon_panels/1_JBLAB_AAprimers_48_basic_panel/',
	'6',	'/scratcha/jblab/amplicon_panels/6_HR_Elkes_panel_2015_08_28/',
	'10',	'/scratcha/jblab/amplicon_panels/10_JBLAB_AAprimers_TP53short_only/',
	'28',	'/scratcha/jblab/amplicon_panels/28_JBLAB_AAprimers_dream_panel/' 
	)

process_amplicon_data = function(amplicon_panel_number) {

	amplicon_panel_path = amplicon_panel_paths %>% 
		dplyr::filter(amplicon_panel_id==amplicon_panel_number) %>%
		dplyr::pull(amplicon_panel_path)

	amplicons_path = paste(amplicon_panel_path,'amplicons.txt', sep='') 
	targets_path = paste(amplicon_panel_path,'targets.txt', sep='') 

	amplicons = readr::read_tsv(
		amplicons_path, 
		comment='@',
		col_names=c('chromosome','amplicon_start','amplicon_stop','strand','amplicon_id')
	)

	targets = readr::read_tsv(
		targets_path, 
		comment='@',
		col_names=c('chromosome','target_start','target_stop','strand','amplicon_id')
	)

	amplicons = amplicons %>% dplyr::select(-strand)
	targets = targets %>% dplyr::select(-strand)

	amplicons = 
		dplyr::inner_join(
			amplicons,
			targets,
			by=c('amplicon_id','chromosome')
		)

	amplicons$amplicon_panel_id = amplicon_panel_number

	return(amplicons)
}

print(amplicon_panel_paths)

all_amplicon_data = purrr::map_dfr(.x=amplicon_panel_paths$amplicon_panel_id, .f=process_amplicon_data)
all_amplicon_data$gene_id = all_amplicon_data$amplicon_id %>% 
	stringr::str_remove('^EXP0116_') %>%
	stringr::str_extract('[A-Z0-9]+') %>%
	stringr::str_remove('EX4$')

all_amplicon_data = all_amplicon_data %>% 
	dplyr::select(amplicon_id, amplicon_panel_id, gene_id, chromosome, amplicon_start, target_start, target_stop, amplicon_stop)

readr::write_tsv(all_amplicon_data, 'amplicon_data_SI.tsv')
