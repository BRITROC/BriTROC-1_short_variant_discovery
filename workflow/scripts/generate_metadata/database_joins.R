# use the database to identify germline samples which were sequenced using amplicon panel 28

library(magrittr)
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

somatic_metadata = dplyr::inner_join(libraries, samples, by=c('fk_sample'='name')) %>%
	dplyr::inner_join(experiments, by=c('fk_experiment'='name')) %>%
	dplyr::inner_join(run_slx, by='fk_slx') %>%
	dplyr::filter(type %in% c('archival','relapse')) %>%
	dplyr::filter(fk_amplicon_panel %in% c(28)) %>%
	dplyr::arrange(fk_britroc_number)

somatic_metadata = somatic_metadata %>% dplyr::mutate(flowcell=stringr::str_extract(string=fk_run, pattern='[-A-Z0-9]+$'))

tp53_somatic_metadata = readr::read_tsv(snakemake@input[['somatic_tp53_metadata']])

patients_with_samples_sequenced_for_tp53 = tp53_somatic_metadata %>% dplyr::pull(fk_britroc_number)

# filter so we only examine non-TP53 genes for patient which had TP53 sequencing for at least one sample
somatic_metadata = somatic_metadata %>% dplyr::filter(fk_britroc_number %in% patients_with_samples_sequenced_for_tp53)

write.table(somatic_metadata, snakemake@output[[1]], row.names=FALSE, quote=FALSE, sep='\t')
