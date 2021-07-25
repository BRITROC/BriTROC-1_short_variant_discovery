
library(magrittr)

sample_level_table = readr::read_tsv('foo.tsv')
patient_level_table = readr::read_tsv('foo2.tsv')

patient_level_table = patient_level_table %>% dplyr::filter(classification=='clonal')
sample_level_table = sample_level_table %>% dplyr::select(-sample_id,-type) %>% unique()

patient_level_table = dplyr::left_join(patient_level_table,sample_level_table, by=c('fk_britroc_number','#Uploaded_variation'))

patient_level_table$Consequence %>% table() %>% print()

print(patient_level_table)
