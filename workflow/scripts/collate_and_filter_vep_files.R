# A simple script within a larger snakemake workflow of collating and filtering variant calls from strelka

library(tidyverse)
library(DBI)
library(RPostgres)

# ensure that the script reads from the users .Renviron text file
readRenviron('~/.Renviron')

vep_files = snakemake@input %>% unlist
sample_names = stringr::str_extract(string=vep_files, pattern='(JBLAB-[0-9]+|IM_[0-9]+)') 

# a rudimentary helper function to add the sample IDs to the annotation output table
mutate_x_y = function(x,y) {
  return(dplyr::mutate(x, sample_id=y))
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
      `#Uploaded_variation` = 'c',
      Location = 'c',
      cDNA_position='c',
      CDS_position='c',
      Protein_position='c',
      Allele='c'
    )
)  %>%
  purrr::map2_dfr(.y=sample_names, .f=mutate_x_y ) %>%
  dplyr::arrange(Gene) %>%
  dplyr::select(sample_id, everything())

# separate columns
extract_column = function(column_name) {

  annotations[[column_name]] = stringr::str_extract(
    string= annotations$Extra,
    pattern = stringr::str_interp('${column_name}=[^;]+'
    )
  )

  annotations[[column_name]] = stringr::str_remove(
    string = annotations[[column_name]],
    pattern = stringr::str_interp('^${column_name}=')
  )
  
  return(
    list(column_name=annotations[[column_name]])
  )
}

new_column_names = c(
		     'IMPACT','SYMBOL','SYMBOL_SOURCE','BIOTYPE','CLIN_SIG','SIFT','PolyPhen','CANONICAL',
                     'ENSP','UNIPARC','EXON','INTRON','SWISSPROT','TREMBL','DOMAINS','FLAGS','gnomAD_AF','AF'
			)

new_columns = purrr::map_dfc(.x=new_column_names, .f=extract_column)
colnames(new_columns) = new_column_names

annotations = cbind(annotations, new_columns) %>% as_tibble()

## begin to filter mutations for predicted pathogenicity

# filter columns
annotations = annotations %>% dplyr::filter(
  !BIOTYPE %in% c('non_stop_decay','nonsense_mediated_decay','processed_pseudogene')
)

annotations = annotations %>% dplyr::filter(
  is.na(gnomAD_AF)
)

annotations = annotations %>% dplyr::filter(
  is.na(AF)
)

annotations = annotations %>% dplyr::filter(
  !grepl('benign',CLIN_SIG)
)

annotations = annotations %>% dplyr::filter(
  !grepl('tolerated',SIFT)
)

annotations = annotations %>% dplyr::filter(
  !grepl('benign',PolyPhen)
)

annotations = annotations %>% dplyr::filter(
  !( grepl('possibly',PolyPhen) & grepl('low_confidence',SIFT))
)

#grepl('possibly',annotations$PolyPhen)
#grepl('low_confidence',annotations$SIFT) 
#c(grepl('possibly',annotations$PolyPhen)) & c(grepl('low_confidence',annotations$SIFT))

#annotations = annotations %>% dplyr::filter(
#  SYMBOL %in% c('BARD1','BRCA1','BRCA2','BRIP1','FANCM','PALB2','RAD51B','RAD51C','RAD51D')
#)

annotations = annotations %>% dplyr::filter(
  CANONICAL == 'YES'
)

annotations = annotations %>% dplyr::filter(
  IMPACT != 'LOW'
)

annotations = annotations %>% dplyr::filter(
  !Consequence %in% c(
    '3_prime_UTR_variant','5_prime_UTR_variant','downstream_gene_variant','intron_variant',
    'upstream_gene_variant', 'non_coding_transcript_exon_variant',
    'intron_variant,non_coding_transcript_variant'
  )
)

# remove transcript isoform duplicates

# this doesn't make a change if we already filter for 'CANONICAL=TRUE'

annotations = annotations %>%
  dplyr::group_by(sample_id,Location,Allele) %>%
  dplyr::filter(dplyr::row_number()==1) %>% dplyr::ungroup()

print(annotations)

annotations %>% group_by(SYMBOL,`#Uploaded_variation`) %>% summarise(n=n()) %>% arrange(n) %>% print()
annotations %>% group_by(SYMBOL) %>% summarise(n=n()) %>% arrange(n) %>% print()

britroc_con <- dbConnect(RPostgres::Postgres(),
                         dbname='britroc1',
                         host='jblab-db.cri.camres.org',
                         port = 5432,
                         user = Sys.getenv('jblab_db_username'),
                         password = Sys.getenv('jblab_db_password')
)

samples = dbReadTable(britroc_con, 'sample')

annotations = dplyr::inner_join(annotations,samples, by=c('sample_id'='name'))

annotations %>% group_by(fk_britroc_number,SYMBOL) %>% summarise(n=n()) %>% ungroup %>% group_by(SYMBOL) %>% summarise(n=n()) %>% arrange(n) %>% print(n=Inf)

write_tsv(annotations, snakemake@output[[1]])

