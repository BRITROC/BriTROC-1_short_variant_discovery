# use the database to identify somatic samples with matched normal samples which were prepared using amplicon panel 28

library(magrittr)
library(DBI)
library(RPostgres)

readRenviron('~/.Renviron')
britroc_con <- DBI::dbConnect(RPostgres::Postgres(),
                 dbname='britroc1',
                 host='jblab-db.cri.camres.org',
                 port = 5432,
                 user = Sys.getenv('jblab_db_username'),
                 password = Sys.getenv('jblab_db_password')
)

clarity_con <- DBI::dbConnect(RPostgres::Postgres(),
                 dbname='clarity',
                 host='jblab-db.cri.camres.org',
                 port = 5432,
                 user = Sys.getenv('jblab_db_username'),
                 password = Sys.getenv('jblab_db_password')
)


source('functions.R')

# example usage
samples = RPostgres::dbReadTable(britroc_con, 'sample')
libraries = RPostgres::dbReadTable(britroc_con, 'slx_library')
experiments = RPostgres::dbReadTable(britroc_con, 'experiment')
run_slx = RPostgres::dbReadTable(britroc_con, 'run_slx')

somatic_metadata = dplyr::inner_join(libraries, samples, by=c('fk_sample'='name')) %>%
	dplyr::inner_join(experiments, by=c('fk_experiment'='name')) %>%
	dplyr::inner_join(run_slx, by='fk_slx', relationship='many-to-many') %>% # many-to-many as some pooled libraries sequenced multiple times
	dplyr::filter(type %in% c('archival','relapse')) %>%
	dplyr::filter(fk_amplicon_panel %in% c(1,10,18,28)) %>%
	dplyr::arrange(fk_britroc_number,type)

somatic_metadata$fk_amplicon_panel = ifelse(somatic_metadata$fk_amplicon_panel == 18, 10, somatic_metadata$fk_amplicon_panel)

somatic_metadata = somatic_metadata %>% dplyr::mutate(flowcell=stringr::str_extract(string=fk_run, pattern='[-A-Z0-9]+$'))

# remove samples which have not passed QC
relevant_samples = remove_non_relevant_samples( 
	non_hgsoc_samples = snakemake@config[['non_hgsoc_samples']],
       samples_with_no_good_sequencing = snakemake@config[['samples_with_no_good_sequencing']],
        samples_with_very_low_purity = snakemake@config[['samples_with_very_low_purity']],
	britroc_con=britroc_con,
	clarity_con=clarity_con,
        analysis_type='cohort'
) %>% dplyr::pull(name)

somatic_metadata = somatic_metadata %>% dplyr::filter(fk_sample %in% relevant_samples)

write.table(somatic_metadata, snakemake@output[[1]], row.names=FALSE, quote=FALSE, sep='\t')
