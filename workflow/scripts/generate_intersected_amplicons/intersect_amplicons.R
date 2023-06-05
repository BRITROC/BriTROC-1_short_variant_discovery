# use the database to identify germline samples which were sequenced using amplicon panel 28

library(magrittr)

panel_6 = readr::read_tsv(snakemake@input[['panel_6']], comment='@', col_names=c('chromosome','start','stop','strand','amplicon'))
panel_28 = readr::read_tsv(snakemake@input[['panel_28']], comment='@', col_names=c('chromosome','start','stop','strand','amplicon'))

# it is important that it is this way round as panel 6 is slightly more restrictive (by one or two bases) for the same amplicon; intervals must be the same for normal and tumour
intersected_amplicons = panel_6 %>% dplyr::filter(amplicon %in% panel_28$amplicon)

print(intersected_amplicons)

write.table(intersected_amplicons, snakemake@output[[1]], row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
