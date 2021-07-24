library(magrittr)

x = readr::read_tsv('results/tumour_sample_vcfs_octopus/124.filtered2.vcf', comment='##')

x = x %>% dplyr::filter(FILTER=='PASS') %>% dplyr::filter(REF=='C') %>% dplyr::filter(ALT=='T')

print(x, width=Inf)
