library(tidyverse)
library(DBI)
library(RPostgres)

# ensure that the script reads from the users .Renviron text file
readRenviron('~/.Renviron')

vep_files = snakemake@input %>% unlist

print(vep_files %>% length)

patient_names = stringr::str_extract(string=vep_files, pattern='[0-9]+') 

print(patient_names)

# a rudimentary helper function to add the sample IDs to the annotation output table
mutate_x_y = function(x,y) {
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
	`#CHROM`='c',
	QUAL='d'
    )
)  %>%
  purrr::map2_dfr(.y=patient_names, .f=mutate_x_y ) %>%
  dplyr::arrange(patient_id) %>%
  dplyr::select(patient_id, everything())

annotations = annotations %>% unique()

print(annotations)



#x = purrr::map_dfr(.x = snakemake@input, .f=readr::read_tsv)

#	library_MAFs = readr::read_tsv(
#                stringr::str_interp('results/tumour_sample_vcfs_octopus/${patient_id}.library_MAFs.vcf')
#                )
#}
