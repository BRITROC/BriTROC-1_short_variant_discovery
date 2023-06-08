# use the database to identify germline samples which were sequenced using amplicon panel 28

library(magrittr)
library(DBI)
library(RPostgres)

readRenviron('~/.Renviron')

# example usage
samples = readr::read_tsv(snakemake@input['DNA_samples_file_path'])
libraries = readr::read_tsv(snakemake@input['slx_library_file_path'])
experiments = readr::read_tsv(snakemake@input['experiments_file_path'])
run_slx = readr::read_tsv(snakemake@input['run_slx_file_path'])

somatic_metadata = dplyr::inner_join(libraries, samples, by=c('fk_sample'='name')) %>%
	dplyr::inner_join(experiments, by=c('fk_experiment'='name')) %>%
	dplyr::inner_join(run_slx, by='fk_slx', relationship='many-to-many') %>%
	dplyr::filter(type %in% c('archival','relapse')) %>%
	dplyr::filter(fk_amplicon_panel %in% c(28)) %>%
	dplyr::arrange(fk_britroc_number)

somatic_metadata = somatic_metadata %>% dplyr::mutate(flowcell=stringr::str_extract(string=fk_run, pattern='[-A-Z0-9]+$'))

tp53_somatic_metadata = readr::read_tsv(snakemake@input['TP53_amplicon_metadata'])

patients_with_samples_sequenced_for_tp53 = tp53_somatic_metadata %>% dplyr::pull(fk_britroc_number)

# filter so we only examine non-TP53 genes for patient which had TP53 sequencing for at least one sample
somatic_metadata = somatic_metadata %>% dplyr::filter(fk_britroc_number %in% patients_with_samples_sequenced_for_tp53)

write.table(somatic_metadata, snakemake@output[[1]], row.names=FALSE, quote=FALSE, sep='\t')
