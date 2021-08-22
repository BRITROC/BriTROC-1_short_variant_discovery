# use the database to identify germline samples which were sequenced using amplicon panel 28
library(magrittr)
library(DBI)
library(RPostgres)

panel_6 = readr::read_tsv(snakemake@input[['panel_6']], comment='@', col_names=c('chromosome','start','stop','strand','amplicon'))
panel_28 = readr::read_tsv(snakemake@input[['panel_28']], comment='@', col_names=c('chromosome','start','stop','strand','amplicon'))

antijoined_amplicons = panel_28 %>% dplyr::filter(!amplicon %in% panel_6$amplicon)

# derive gene symbols from amplicon identifiers
antijoined_amplicons = antijoined_amplicons %>% dplyr::mutate(gene_symbol = amplicon %>% stringr::str_extract('[A-Z]+[0-9]*[A-Z]*')) 

# remove gene symbols which are found in both panel 28 and panel 6
antijoined_amplicons = antijoined_amplicons %>% dplyr::filter(!gene_symbol %in% c('TP53','BRCA2','EXP0116'))
antijoined_amplicons = antijoined_amplicons %>% dplyr::select(-gene_symbol)

# order amplicons by chromosome
antijoined_amplicons = antijoined_amplicons %>% 
	dplyr::mutate(chromosome_integer = chromosome %>% stringr::str_remove('chr') %>% as.integer) %>% 
	dplyr::arrange(chromosome_integer) %>%
	dplyr::select(-chromosome_integer)

write.table(antijoined_amplicons, snakemake@output[[1]], row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
