# use the database to identify somatic samples with matched normal samples which were prepared using amplicon panel 28

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
	dplyr::filter(fk_amplicon_panel %in% c(28)) %>%
	dplyr::arrange(fk_britroc_number,type)

germline_metadata = germline_metadata %>% mutate(flowcell=stringr::str_extract(string=fk_run, pattern='[-A-Z0-9]+$'))

somatic_metadata = dplyr::inner_join(libraries, samples, by=c('fk_sample'='name')) %>%
	dplyr::inner_join(experiments, by=c('fk_experiment'='name')) %>%
	dplyr::inner_join(run_slx, by='fk_slx') %>%
	dplyr::filter(type %in% c('archival','relapse')) %>%
	dplyr::filter(fk_britroc_number %in% germline_metadata$fk_britroc_number) %>%
	dplyr::filter(fk_amplicon_panel %in% c(28)) %>%
	dplyr::arrange(fk_britroc_number,type)

somatic_metadata = somatic_metadata %>% mutate(flowcell=stringr::str_extract(string=fk_run, pattern='[-A-Z0-9]+$'))

write.table(somatic_metadata, snakemake@output[['somatic_metadata']], row.names=FALSE, quote=FALSE, sep='\t')
