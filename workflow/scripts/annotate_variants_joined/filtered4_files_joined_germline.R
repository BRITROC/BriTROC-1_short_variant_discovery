# A simple script within a larger snakemake workflow of collating and filtering variant calls outputted by a variant calling algorithm

library(tidyverse)
library(DBI)
library(RPostgres)

# ensure that the script reads from the users .Renviron text file
readRenviron('~/.Renviron')

vep_files = snakemake@input %>% unlist

patient_names = stringr::str_extract(string=vep_files, pattern='[0-9]+.germline.filtered.vcf$') %>%
        stringr::str_extract('[0-9]+')

# a rudimentary helper function to add the sample IDs to the annotation output table
mutate_x_y = function(x,y) {
  
  x = x %>% dplyr::select(-dplyr::matches('JBLAB-'))
  x = x %>% dplyr::select(-dplyr::matches('IM_'))

  return(dplyr::mutate(x, patient_id=y))
}

# A pipe which reads in each VEP file, adds the sample ID as an additional column and reformats
annotations = purrr::map(
  .x=vep_files,
  .f=readr::read_tsv, 
  comment='##', 
  skip=0, 
  guess_max = 5000,
  na = '-',
  col_names=TRUE,
  col_types =
    readr::cols(
	POS='i',
	REF='c',
	ALT='c',
	`#CHROM`='c'
    )
)  %>%
  purrr::map2_dfr(.y=patient_names, .f=mutate_x_y ) %>%
  dplyr::arrange(patient_id) %>%
  dplyr::select(patient_id, everything())

annotations = annotations %>% unique()

readr::write_tsv(annotations, path=snakemake@output[[1]], append=FALSE)
quit()
