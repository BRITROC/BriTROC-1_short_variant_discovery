# use the database to identify germline samples which were sequenced using amplicon panel 28

library(magrittr)
library(DBI)
library(RPostgres)

readRenviron('~/.Renviron')

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

write.table(germline_metadata, snakemake@output[[1]], row.names=FALSE, quote=FALSE, sep='\t')
