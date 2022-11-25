# use the database to identify germline samples which were sequenced using amplicon panel 28

library(tidyverse)
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

# example usage
samples = RPostgres::dbReadTable(britroc_con, 'sample')
libraries = RPostgres::dbReadTable(britroc_con, 'slx_library')
experiments = RPostgres::dbReadTable(britroc_con, 'experiment')
run_slx = RPostgres::dbReadTable(britroc_con, 'run_slx')

germline_metadata = dplyr::inner_join(libraries, samples, by=c('fk_sample'='name')) %>%
	dplyr::inner_join(experiments, by=c('fk_experiment'='name')) %>%
	dplyr::inner_join(run_slx, by='fk_slx') %>%
	dplyr::filter(type %in% c('germline')) %>%
	dplyr::filter(fk_amplicon_panel %in% c(6,28)) %>%
	dplyr::arrange(fk_britroc_number)

germline_metadata = germline_metadata %>% mutate(flowcell=stringr::str_extract(string=fk_run, pattern='[-A-Z0-9]+$'))

somatic_metadata = dplyr::inner_join(libraries, samples, by=c('fk_sample'='name')) %>%
	dplyr::inner_join(experiments, by=c('fk_experiment'='name')) %>%
	dplyr::inner_join(run_slx, by='fk_slx') %>%
	dplyr::filter(type %in% c('archival','relapse')) %>%
	dplyr::filter(fk_britroc_number %in% germline_metadata$fk_britroc_number) %>% # ensure there is a relevant germline sample
	dplyr::filter(fk_amplicon_panel %in% c(6,28)) %>%
	dplyr::arrange(fk_britroc_number)

somatic_metadata = somatic_metadata %>% mutate(flowcell=stringr::str_extract(string=fk_run, pattern='[-A-Z0-9]+$'))

somatic_metadata %>% tibble::as_tibble() %>% print()

somatic_metadata = somatic_metadata %>% dplyr::filter(fk_britroc_number %in% germline_metadata$fk_britroc_number)
germline_metadata = germline_metadata %>% dplyr::filter(fk_britroc_number %in% somatic_metadata$fk_britroc_number)

write.table(germline_metadata, snakemake@output[['germline_metadata']], row.names=FALSE, quote=FALSE, sep='\t')
write.table(somatic_metadata, snakemake@output[['somatic_metadata']], row.names=FALSE, quote=FALSE, sep='\t')