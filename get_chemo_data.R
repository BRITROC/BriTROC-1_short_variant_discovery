source('~/.Renviron')
source('functions.R')

library(DBI)
library(RPostgres)
library(magrittr)

britroc_con = make_connection_to_postgres_server('britroc1', 'jblab-db.cri.camres.org', 5432)

chemo_table = dbReadTable(britroc_con, 'chemotherapy_lines') %>% tibble::as_tibble()
chemo_table = chemo_table %>% dplyr::filter(time_relation_to_image_guided_biopsy=='before')

chemo_drugs_table = dbReadTable(britroc_con, 'chemotherapy_lines_drugs') %>% tibble::as_tibble()

chemo_drugs_table = chemo_drugs_table %>%
	dplyr::semi_join(
	chemo_table,
	by=c('fk_britroc_number','chemotherapy_line')
	)

print(chemo_table)
print(chemo_drugs_table)

chemo_drugs_table = chemo_drugs_table %>% dplyr::select(-id)

readr::write_tsv(chemo_table, 'BriTROC-1_chemotherapy_lines_before_study_entry.tsv')
readr::write_tsv(chemo_drugs_table, 'BriTROC-1_chemotherapy_lines_drugs_before_study_entry.tsv')
