# use the database to identify germline samples which were sequenced using amplicon panel 28

library(tidyverse)
library(DBI)
library(RPostgres)

panel_6 = readr::read_tsv(snakemake@input[['panel_6']], comment='@', col_names=c('chromosome','start','stop','strand','amplicon'))
panel_28 = readr::read_tsv(snakemake@input[['panel_28']], comment='@', col_names=c('chromosome','start','stop','strand','amplicon'))

antijoined_amplicons = panel_28 %>% filter(!amplicon %in% panel_6$amplicon)

write.table(antijoined_amplicons, snakemake@output[[1]], row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
