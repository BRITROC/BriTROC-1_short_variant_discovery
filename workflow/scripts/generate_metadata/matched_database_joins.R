# use the database to identify germline samples which were sequenced using amplicon panel 28

library(magrittr)
library(DBI)
library(RPostgres)

samples = readr::read_tsv(snakemake@input['DNA_samples_file_path'])
libraries = readr::read_tsv(snakemake@input['slx_library_file_path'])
experiments = readr::read_tsv(snakemake@input['experiments_file_path'])
run_slx = readr::read_tsv(snakemake@input['run_slx_file_path'])

germline_metadata = dplyr::inner_join(libraries, samples, by=c('fk_sample'='name')) %>%
	dplyr::inner_join(experiments, by=c('fk_experiment'='name')) %>%
	dplyr::inner_join(run_slx, by='fk_slx', relationship='many-to-many') %>%
	dplyr::filter(type %in% c('germline')) %>%
	dplyr::filter(fk_amplicon_panel %in% c(6,28)) %>%
	dplyr::arrange(fk_britroc_number)

germline_metadata = germline_metadata %>% dplyr::mutate(flowcell=stringr::str_extract(string=fk_run, pattern='[-A-Z0-9]+$'))

somatic_metadata = dplyr::inner_join(libraries, samples, by=c('fk_sample'='name')) %>%
	dplyr::inner_join(experiments, by=c('fk_experiment'='name')) %>%
	dplyr::inner_join(run_slx, by='fk_slx', relationship='many-to-many') %>%
	dplyr::filter(type %in% c('archival','relapse')) %>%
	dplyr::filter(fk_britroc_number %in% germline_metadata$fk_britroc_number) %>% # ensure there is a relevant germline sample
	dplyr::filter(fk_amplicon_panel %in% c(6,28)) %>%
	dplyr::arrange(fk_britroc_number)

somatic_metadata = somatic_metadata %>% dplyr::mutate(flowcell=stringr::str_extract(string=fk_run, pattern='[-A-Z0-9]+$'))

# ensure that there is a archival and a relapse sample for the same patient
paired_somatic_metadata = somatic_metadata %>% 
	dplyr::group_by(fk_britroc_number,type) %>% 
	dplyr::summarise(n=dplyr::n()) %>% 
	dplyr::ungroup() %>% 
	dplyr::group_by(fk_britroc_number) %>% 
	dplyr::summarise(n=dplyr::n()) %>%
	dplyr::filter(n>1) %>%
	dplyr::select(-n)

somatic_metadata = somatic_metadata %>% dplyr::filter(fk_britroc_number %in% paired_somatic_metadata$fk_britroc_number)

somatic_metadata = somatic_metadata %>% dplyr::filter(fk_britroc_number %in% germline_metadata$fk_britroc_number)
somatic_metadata = somatic_metadata %>% dplyr::select(fk_britroc_number) %>% unique()

write.table(somatic_metadata, snakemake@output[[1]], row.names=FALSE, quote=FALSE, sep='\t')
