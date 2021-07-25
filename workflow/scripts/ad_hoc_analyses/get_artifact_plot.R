library(ggplot2)
library(magrittr)
library(DBI)
library(RPostgres)

num_artifacts = purrr::map_dfr(.x=snakemake@input, readr::read_tsv)

metadata = readr::read_tsv('config/matched_somatic_metadata.tsv')

relapse_data = readr::read_csv('resources/Biopsies_relapse.csv') %>% dplyr::select(JBLAB_ID, 'NBF/UMFIX', 'DNA - Qubit [ng/ul]','Total DNA (ug)')
relapse_data = relapse_data %>% dplyr::rename('fixation_method'='NBF/UMFIX')
relapse_data$`DNA - Qubit [ng/ul]` = as.double(relapse_data$`DNA - Qubit [ng/ul]`)
relapse_data$`Total DNA (ug)` = as.double(relapse_data$`Total DNA (ug)`)

print(relapse_data)

num_artifacts = num_artifacts %>% dplyr::inner_join(metadata, by=c('sample_id'='fk_sample')) %>% dplyr::select(sample_id,num_artifacts,MAF_upper_quartile) %>% unique()

readRenviron('~/.Renviron')

britroc_con <- dbConnect(RPostgres::Postgres(),
                 dbname='britroc1',
                 host='jblab-db.cri.camres.org', # use 'jblab-db.cruk.cam.ac.uk' instead for external IPs,
                 port = 5432,
                 user = Sys.getenv('jblab_db_username'),
                 password = Sys.getenv('jblab_db_password')
)

patients = dbReadTable(britroc_con, 'patients')
samples = dbReadTable(britroc_con, 'sample')

print(samples %>% tibble::as_tibble())

num_artifacts = num_artifacts %>% dplyr::inner_join(samples, by=c('sample_id'='name')) %>% 
	dplyr::select(sample_id,num_artifacts,fk_britroc_number,type,MAF_upper_quartile) %>% unique()

num_artifacts = num_artifacts %>% dplyr::inner_join(patients, by=c('fk_britroc_number'='britroc_number')) %>% 
	dplyr::select(sample_id,num_artifacts,type,fk_britroc_number,site,MAF_upper_quartile) %>% unique()

num_artifacts = num_artifacts %>% 
	dplyr::left_join(relapse_data, by=c('sample_id'='JBLAB_ID'))

num_artifacts$`fixation_method` = ifelse(num_artifacts$type == 'archival', 'archival', num_artifacts$`fixation_method`)
num_artifacts$`fixation_method` = num_artifacts$`fixation_method` %>% tidyr::replace_na(replace='unknown')

ggplot(num_artifacts, aes(x=site, y=num_artifacts)) + 
	geom_violin() + 
	geom_point(position = position_jitter(width=0.20), aes(colour=`fixation_method`)) +
	facet_wrap(vars(type), dir='v') + ylab('Number of artifacts')

print(num_artifacts, width=Inf)

ggsave('artifact_by_site.png', device='png', width=32)

ggplot(num_artifacts, aes(x=site, y=MAF_upper_quartile)) + 
	geom_violin() + 
	geom_point(position = position_jitter(width=0.20), aes(colour=`fixation_method`)) +
	facet_wrap(vars(type), dir='v') + ylab('95th quantile of artifact MAF distribution')

ggsave('artifact_by_site2.png', device='png', width=32)

#ggplot(num_artifacts, aes(x=`DNA - Qubit [ng/ul]`, y=num_artifacts)) + 
#	geom_point(aes(colour=`fixation_method`)) +
#	facet_wrap(vars(type), dir='v')

#ggsave('artifacts_by_conc.png', device='png', width=32)

#ggplot(num_artifacts, aes(x=`Total DNA (ug)`, y=num_artifacts)) + 
#	geom_point(aes(colour=`fixation_method`)) +
#	facet_wrap(vars(type), dir='v')

#ggsave('artifacts_by_total_dna.png', device='png', width=32)

#num_artifacts = num_artifacts %>% 
#	dplyr::left_join(relapse_data, by=c('sample_id'='JBLAB_ID')) %>% 
#	dplyr::filter(! (is.na(`fixation_method`) & type=='relapse') )
#
#num_artifacts$`fixation_method` = dplyr::coalesce(num_artifacts$`fixation_method`, num_artifacts$type)

print(num_artifacts %>% dplyr::filter(`fixation_method`=='NBF'), n=Inf)

print(num_artifacts %>% dplyr::filter(grepl('Bartholomew',site), type=='relapse') %>% .$`fixation_method` %>% table(), n=Inf)
print(num_artifacts %>% dplyr::filter(grepl('Beatson',site), type=='relapse') %>% .$`fixation_method` %>% table(), n=Inf)
print(num_artifacts %>% dplyr::filter(grepl('Hammersmith',site), type=='relapse') %>% .$`fixation_method` %>% table(), n=Inf)

num_artifacts$`fixation_method` %>% table()

ggplot(num_artifacts, aes(x=`fixation_method`, y=num_artifacts)) + geom_violin() + geom_point(position = position_jitter(width=0.20)) +
	ylab('number of artifacts')

ggsave(snakemake@output[[1]], device='png')

ggplot(num_artifacts, aes(x=`fixation_method`, y=MAF_upper_quartile)) + geom_violin() + geom_point(position = position_jitter(width=0.20)) +
	ylab('95th quantile of artifact MAF distribution')

ggsave('article_by_fixation_method.png', device='png')
