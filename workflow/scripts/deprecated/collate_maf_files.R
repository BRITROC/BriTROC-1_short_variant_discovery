# A simple script used to collate maf files together

#R.Version()
#.libPaths()
library(tidyverse)

sample_ids = stringr::str_extract(
  string=snakemake@input,
  pattern='[^\\.]+')

annotations = purrr::map(
  .x=snakemake@input,
  .f=read_tsv, comment='##', skip=1, guess_max = 5000,
  col_types = cols(Transcript_Position = 'i', RU = 'c', Chromosome='c', HGNC_Gene_Family_ID='c')
) %>% map2(.y=sample_ids, .f=function(x,y) (mutate(x, sample_id=y)) ) %>%
  plyr::ldply(data.frame, check.name=FALSE) %>% arrange(Hugo_Symbol) %>% select(sample_id, everything())

print(annotations)

#purrr::map(.x=snakemake@input, .f=readr::read_tsv, comment='##', skip=1) %>% purrr::map(.f=map_dfr)
