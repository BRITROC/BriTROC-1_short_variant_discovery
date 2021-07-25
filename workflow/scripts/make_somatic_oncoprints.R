# A scrip to generate oncoprints from somatic short variant data

library(tidyverse)
library(DBI)
library(RPostgres)
library(ComplexHeatmap)

archival_variants = readr::read_tsv('filtered_archival_vep_calls_octopus.tsv')
archival_variants_tp53 = readr::read_tsv('filtered_archival_vep_calls_octopus_tp53.tsv')

archival_variants = rbind(archival_variants, archival_variants_tp53)

archival_variants$Consequence = factor(archival_variants$Consequence, 
                                       levels=c(
                                         'frameshift_variant',
                                         'inframe_deletion',
                                         'inframe_insertion',
                                         'stop_gained',
                                         'stop_gained,frameshift_variant',
                                         'stop_gained,splice_region_variant',
                                         'missense_variant',
                                         'splice_donor_variant',
                                         'splice_acceptor_variant',
                                         'splice_region_variant,intron_variant'
                                         ))

archival_variants = archival_variants %>% arrange(fk_britroc_number,Consequence)

archival_variants = archival_variants %>% filter(`#Uploaded_variation`!='chr17_41231352_T/C')

# retrieve all relevant samples

britroc_con <- dbConnect(RPostgres::Postgres(),
                         dbname='britroc1',
                         host='jblab-db.cri.camres.org',
                         port = 5432,
                         user = Sys.getenv('jblab_db_username'),
                         password = Sys.getenv('jblab_db_password')
)

clarity_con <- dbConnect(RPostgres::Postgres(),
                         dbname='clarity',
                         host='jblab-db.cri.camres.org',
                         port = 5432,
                         user = Sys.getenv('jblab_db_username'),
                         password = Sys.getenv('jblab_db_password')
)

samples = dbReadTable(britroc_con, 'sample')

slx_library = dbReadTable(britroc_con, 'slx_library') %>%
       dplyr::filter(fk_slx!='SLX-13716') %>% # no data for this SLX
       dplyr::filter(grepl('AA',fk_experiment)) # select for tam-seq experiments only

slx_clarity = dbReadTable(clarity_con, 'slx')

slx_library = dplyr::semi_join(slx_library, slx_clarity, by=c('fk_slx'='name')) # ensure slx is actually in clarity

experiments = dbReadTable(britroc_con, 'experiment') %>% dplyr::filter(fk_amplicon_panel %in% c(6,28)) # only allow panel 6 or panel 28 amplicon panels
slx_library = slx_library %>% dplyr::semi_join(experiments, by=c('fk_experiment'='name'))

# retrieve analysis type
#analysis_type = stringr::str_extract(string=snakemake@output[[1]], pattern='(archival|relapse)')
#print(analysis_type)

sequenced_samples = dplyr::semi_join(samples, slx_library, by=c('name'='fk_sample'))
sequenced_samples = sequenced_samples %>%
       dplyr::group_by(fk_britroc_number,type) %>%
       dplyr::summarise(n=n()) %>%
       tidyr::pivot_wider(names_from=type, values_from=n) %>%
       dplyr::filter(!is.na(archival) & !is.na(relapse) & !is.na(germline)) # https://stackoverflow.com/questions/27197617/filter-data-frame-by-character-column-name-in-dplyr

#num_patients_with_paired_samples = sequenced_samples %>% .$fk_britroc_number %>% as.factor %>% levels %>% length()

archival_variants = 
  archival_variants %>% dplyr::select(fk_britroc_number,SYMBOL, `#Uploaded_variation`, Existing_variation, Consequence) %>%
  unique()

## read in TP53 mutations
#tp53_mutations = readr::read_tsv('sample_sheet.tsv')
#tp53_mutations = tp53_mutations %>% dplyr::semi_join(
#  samples %>% filter(type=='archival'), 
#  by=c('SAMPLE_ID'='name'))

#tp53_mutations = tp53_mutations %>% 
#  select(PATIENT_ID,TP53freq) %>% 
#  filter(!is.na(TP53freq)) %>%
#  select(PATIENT_ID) %>%
#  unique() %>%
#  rename(fk_britroc_number=PATIENT_ID) %>%
#  mutate(
#    SYMBOL='TP53', 
#    `#Uploaded_variation`=NA, 
#    Existing_variation=NA, 
#    Consequence='unknown'
#    )

#tp53_mutations$fk_britroc_number = stringr::str_remove(string=tp53_mutations$fk_britroc_number, pattern='BRITROC-')
#tp53_mutations$fk_britroc_number = as.integer(tp53_mutations$fk_britroc_number)

#subset for only those patients we are interested in
#tp53_mutations = tp53_mutations %>% filter(fk_britroc_number %in% sequenced_samples$fk_britroc_number)

#archival_variants = rbind(archival_variants,tp53_mutations)

# remove synonymous and intron variants
archival_variants = 
  archival_variants %>% dplyr::filter(!Consequence %in% c('intron_variant','synonymous_variant')) %>%
  unique()

somatic_samples_with_no_mutations = dplyr::anti_join(
  sequenced_samples, archival_variants, by=c('fk_britroc_number')
)

# reformat
somatic_samples_with_no_mutations = somatic_samples_with_no_mutations %>%
  select(fk_britroc_number) %>% 
  unique() %>% 
  mutate(variant_type='dummy', gene_symbol='BRCA1') %>% 
  select(fk_britroc_number,variant_type,gene_symbol)

somatic_samples_with_mutations = archival_variants %>% select(fk_britroc_number,SYMBOL,Consequence)
somatic_samples_with_no_mutations =  somatic_samples_with_no_mutations %>% ungroup()
colnames(somatic_samples_with_mutations) = c('fk_britroc_number','gene_symbol','variant_type')

somatic_variants = rbind(somatic_samples_with_mutations, somatic_samples_with_no_mutations)

# remap variant type categories
somatic_variants$variant_type = dplyr::recode(
  somatic_variants$variant_type,
  'frameshift,splice_region'='frameshift',
  'stop_gained,frameshift_variant'='frameshift',
  'frameshift_variant'='frameshift',
  'inframe_deletion'='inframe indel',
  'inframe_insertion'='inframe indel',
  'stop_gained,splice_region_variant'='stop_gained',
  'splice_acceptor'='splice_region_SNV',
  'splice_region,intron'='splice_region_SNV',
  'splice_region,synonymous'='splice_region_SNV',
  'splice_acceptor_variant'='splice_region_SNV',
  'splice_donor_variant'='splice_region_SNV',
  'nonsynonymous'='missense',
  'missense_variant'='missense',
  "splice_region_variant,intron_variant" = 'splice_region_SNV',
  "splice_region_variant,synonymous_variant"='splice_region_SNV',
  'splice_region_variant'='splice_region_SNV',
  'dummy'='' # make the dummy row blank
)

# for simplicity for now - only one mutation per patient-gene group
somatic_variants = somatic_variants %>% as_tibble() %>% 
  dplyr::group_by(fk_britroc_number,gene_symbol) %>%
  dplyr::filter(dplyr::row_number()==1) %>% dplyr::ungroup()

# convert the data frame to a matrix
somatic_variants = somatic_variants %>% 
  pivot_wider(names_from=fk_britroc_number, values_from=variant_type) %>% 
  as.matrix()

row.names(somatic_variants) = somatic_variants[,1]
somatic_variants = somatic_variants[,-1]

# replace NA values
somatic_variants = replace_na(somatic_variants, replace = '')

# map short mutation classes to labels
col = c(
  'frameshift'= 'red',
  'stop_gained'= 'orange',
  'splice_region_SNV'= 'yellow',
  'inframe indel'= 'cyan',
  'missense'= 'blue'
)

alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),   
  frameshift = alter_graphic("rect", fill = col["frameshift"]),
  stop_gained = alter_graphic("rect", fill = col["stop_gained"]),
  `inframe indel` = alter_graphic("rect", fill = col["inframe indel"]),
  splice_region_SNV = alter_graphic("rect", fill = col["splice_region_SNV"]),
  missense = alter_graphic("rect", fill = col["missense"])
)

ComplexHeatmap::oncoPrint(
  mat=somatic_variants,
  alter_fun = alter_fun,
  heatmap_legend_param = list(labels=c('frameshift','stop gained','inframe indel','splice region SNV','missense'),
                              title='Alterations')
)
