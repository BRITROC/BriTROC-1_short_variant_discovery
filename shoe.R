
library(magrittr)

archival_variants = readr::read_tsv('tmp_annotations_joined_archival.tsv') %>% dplyr::select(patient_id,CHROM,POS,REF,ALT)

print(archival_variants)

metadata = readr::read_tsv('config/matched_somatic_metadata.tsv') %>% dplyr::filter(type=='archival') %>% dplyr::select(fk_britroc_number,fk_slx)

archival_variants = dplyr::inner_join(archival_variants, metadata, by=c('patient_id'='fk_britroc_number')) %>% unique() %>%
	dplyr::arrange(patient_id,CHROM,fk_slx)

print(archival_variants)
archival_variants$fk_slx %>% table()

readr::write_tsv(archival_variants, 'shoe.txt')
