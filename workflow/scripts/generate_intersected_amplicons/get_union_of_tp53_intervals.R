library(tidyverse)

panel_1 = readr::read_tsv(snakemake@input[['panel_1']], comment='@', col_names = c('chromosome', 'start','end','strand','amplicon'))
panel_10 = readr::read_tsv(snakemake@input[['panel_10']], comment='@', col_names = c('chromosome', 'start','end','strand','amplicon'))
panel_28 = readr::read_tsv(snakemake@input[['panel_28']], comment='@', col_names = c('chromosome', 'start','end','strand','amplicon'))

all_panels = rbind(panel_10,panel_28) %>% unique()
all_panels = all_panels %>% dplyr::filter(grepl('(TP53|tp53)',amplicon))
all_panels = all_panels %>% arrange(start,end,amplicon)

print(all_panels)

write_tsv(all_panels, snakemake@output[[1]], col_names=FALSE)
