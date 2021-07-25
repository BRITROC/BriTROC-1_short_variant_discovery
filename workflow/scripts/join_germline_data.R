library(tidyverse)
library(DBI)
library(RPostgres)

britroc_con <- dbConnect(RPostgres::Postgres(),
                         dbname='britroc1',
                         host='jblab-db.cri.camres.org',
                         port = 5432,
                         user = Sys.getenv('jblab_db_username'),
                         password = Sys.getenv('jblab_db_password')
)

germline_snvs = dbReadTable(britroc_con, 'germline_snvs') %>% 
	as_tibble() %>%
	select(fk_sample,chromosome,position,ref,alt)
#germline_indels = dbReadTable(britroc_con, 'germline_indels') %>% 
#	as_tibble() %>%
#	select(fk_sample,chromosome,position,ref,alt)

germline_variants = rbind(germline_snvs)  #,germline_indels)

samples = dbReadTable(britroc_con, 'sample') %>%
	as_tibble()

germline_variants = inner_join(germline_variants,samples, by=c('fk_sample'='name')) %>% unique()

#print(germline_variants)

archival_variants = readr::read_tsv('results/filtered_archival_vep_calls.tsv') %>% 
	select(sample_id,fk_britroc_number,`#Uploaded_variation`)

#print(archival_variants)

#print(germline_variants$ref %>% as.factor %>% levels())
germline_variants = germline_variants %>% unite(col='join_column_tmp', chromosome,position,ref, sep='_', remove=TRUE) %>%
	unite(col='join_column', join_column_tmp, alt, sep='/', remove=TRUE)


semi_join(archival_variants, germline_variants, by=c('#Uploaded_variation'='join_column'))
	 
germline_variants %>% filter(fk_britroc_number==4)
archival_variants %>% filter(fk_britroc_number==4)
