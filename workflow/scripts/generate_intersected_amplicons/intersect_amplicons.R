# use the database to identify germline samples which were sequenced using amplicon panel 28

library(tidyverse)
library(DBI)
library(RPostgres)

panel_6 = readr::read_tsv(snakemake@input[['panel_6']], comment='@', col_names=c('chromosome','start','stop','strand','amplicon'))
panel_28 = readr::read_tsv(snakemake@input[['panel_28']], comment='@', col_names=c('chromosome','start','stop','strand','amplicon'))

# it is important that it is this way round as panel 6 is slightly more restrictive (by one or two bases) for the same amplicon; intervals must be the same for normal and tumour
intersected_amplicons = panel_6 %>% filter(amplicon %in% panel_28$amplicon)

print(intersected_amplicons)

write.table(intersected_amplicons, snakemake@output[[1]], row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

#britroc_con <- dbConnect(RPostgres::Postgres(),
#                 dbname='britroc1',
#                 host='jblab-db.cri.camres.org',
#                 port = 5432,
#                 user = Sys.getenv('jblab_db_username'),
#                 password = Sys.getenv('jblab_db_password')
#)
#
## read in the data
#amplicons = RPostgres::dbReadTable(britroc_con, 'amplicon')
#amplicon_panel_amplicon = RPostgres::dbReadTable(britroc_con, 'amplicon_panel_amplicon')
#gene = RPostgres::dbReadTable(britroc_con, 'gene')
#
## intersect amplicons
#panel_28 = amplicon_panel_amplicon %>% filter(fk_amplicon_panel==28)
#panel_6 = amplicon_panel_amplicon %>% filter(fk_amplicon_panel==6)
#amplicons = amplicons %>% filter(name %in% panel_28$fk_amplicon) %>% filter(name %in% panel_6$fk_amplicon)  
#
## join gene chromosome information
#amplicons = inner_join(amplicons, gene, by=c('fk_gene'='name')) 
#amplicons$chromosome = paste('chr',amplicons$chromosome, sep='')
#
#amplicons$score = 1
#
## convert to bed format
#amplicons = amplicons %>% select(chromosome,start,stop,name,score,)
#
#print(amplicons %>% as_tibble)
#
## write to disk
##write.table(germline_metadata, snakemake@output[['germline_metadata']], row.names=FALSE, quote=FALSE, sep='\t')
