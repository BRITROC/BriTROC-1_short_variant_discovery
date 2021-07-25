library(magrittr)
library(DBI)
library(RPostgres)
library(ComplexHeatmap)

#archival_non_tp53 = readr::read_tsv('results/filtered_archival_vep_calls_octopus_joined.tsv')
relapse_non_tp53 = readr::read_tsv('results/filtered_relapse_vep_calls_octopus_joined.tsv')
relapse_tp53 = readr::read_tsv('foo.tsv') %>% dplyr::filter(type=='relapse') %>% dplyr::mutate(SYMBOL='TP53')

tp53_variant_clonality = readr::read_tsv('foo2.tsv')

relapse_tp53_clonal_mutations = dplyr::inner_join(
	relapse_tp53, 
	tp53_variant_clonality %>% dplyr::filter(classification=='clonal' | is.na(classification)), 
	by=c('fk_britroc_number','#Uploaded_variation')
) %>% dplyr::select(fk_britroc_number,`#Uploaded_variation`,Consequence,SYMBOL)

relapse_tp53_no_clonal_mutation = relapse_tp53 %>% 
	dplyr::filter(!fk_britroc_number %in% relapse_tp53_clonal_mutations$fk_britroc_number) %>%
	dplyr::select(fk_britroc_number,`#Uploaded_variation`,Consequence,SYMBOL)

relapse_tp53 = rbind(relapse_tp53_clonal_mutations, relapse_tp53_no_clonal_mutation)

#print(relapse_tp53)

relapse_non_tp53 = relapse_non_tp53 %>% dplyr::select(patient_id,`#Uploaded_variation`,Consequence,SYMBOL)
relapse_tp53 = relapse_tp53 %>% dplyr::rename(patient_id=fk_britroc_number)

relapse_non_tp53 = relapse_non_tp53 %>% dplyr::filter(Consequence %in% c('frameshift_variant','stop_gained'))
relapse_non_tp53 = relapse_non_tp53 %>% dplyr::filter(SYMBOL %in% c('BRCA1','BRCA2','FANCM','BARD1'))

print(relapse_non_tp53 %>% dplyr::select(patient_id, `#Uploaded_variation`, SYMBOL, Consequence), n=Inf)

relapse_variants = rbind(relapse_tp53, relapse_non_tp53)

relapse_variants$Consequence = factor(relapse_variants$Consequence,
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
                                         'splice_region_variant,intron_variant',
					 'splice_region_variant,synonymous_variant'
                                         ))

relapse_variants = relapse_variants %>% dplyr::arrange(patient_id,Consequence)

relapse_variants = relapse_variants %>% tibble::add_row(SYMBOL='FANCM', patient_id=176, Consequence='frameshift')

relapse_variants = relapse_variants %>% 
	dplyr::filter(`#Uploaded_variation`!='chr17_41231352_T/C') %>%
	dplyr::filter(!(SYMBOL == 'FANCM' & patient_id==123)) %>%
	dplyr::filter(!(SYMBOL == 'BRCA2' & patient_id==77))

# retrieve sequencing metadata

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
#analysis_type = stringr::str_extract(string=snakemake@output[[1]], pattern='(relapse|relapse)')
#print(analysis_type)

sequenced_samples = dplyr::semi_join(samples, slx_library, by=c('name'='fk_sample'))

non_hgsoc_samples = c('JBLAB-4114','JBLAB-4916','IM_249','IM_250','IM_234','IM_235','IM_236','IM_237','JBLAB-4271','IM_420',
			'IM_262','JBLAB-4922','JBLAB-4923','IM_303','IM_290','IM_43','IM_293','IM_307','IM_308','IM_309','IM_424',
			'IM_302','IM_303','IM_304','IM_305',
			'JBLAB-19320','IM_61','IM_62','IM_63','IM_397','IM_302','IM_98','JBLAB-4210','IM_147','JBLAB-4216','IM_44')

samples_with_no_good_sequencing = c('IM_144','IM_435','IM_436','IM_158','IM_296','IM_373','IM_154','IM_297','IM_365','IM_432','IM_429','IM_368','IM_441') 

samples_with_very_low_purity = c('IM_1','IM_2','IM_3','IM_4','IM_20','IM_26','IM_27','IM_69','IM_86','IM_90','IM_93','IM_94','IM_173','IM_177','IM_179','IM_200',
				'IM_241','IM_242','IM_417','IM_418','IM_419','IM_420','IM_221','IM_264','IM_329','IM_289','IM_308','IM_309',
				'IM_338','IM_339','IM_340','IM_341','IM_342','IM_432','IM_372','IM_272','IM_392')

bad_samples = c(non_hgsoc_samples, samples_with_no_good_sequencing, samples_with_very_low_purity)

sequenced_samples = sequenced_samples %>% dplyr::filter(!name %in% bad_samples)

sequenced_samples = sequenced_samples %>%
       dplyr::group_by(fk_britroc_number,type) %>%
       dplyr::summarise(n=dplyr::n()) %>%
       tidyr::pivot_wider(names_from=type, values_from=n) %>%
       dplyr::filter(!is.na(archival) & !is.na(relapse) & !is.na(germline)) # https://stackoverflow.com/questions/27197617/filter-data-frame-by-character-column-name-in-dplyr

relapse_variants = relapse_variants %>% dplyr::filter(patient_id %in% sequenced_samples$fk_britroc_number)

#####

relapse_variants =
  relapse_variants %>% dplyr::select(patient_id,SYMBOL, `#Uploaded_variation`, Consequence) %>%
  unique()

# remove synonymous and intron variants
relapse_variants = 
  relapse_variants %>% dplyr::filter(!Consequence %in% c('intron_variant','synonymous_variant')) %>%
  unique()

somatic_samples_with_no_mutations = dplyr::anti_join(
  sequenced_samples, relapse_variants, by=c('fk_britroc_number'='patient_id')
)

somatic_samples_with_no_mutations = somatic_samples_with_no_mutations %>% dplyr::rename(patient_id='fk_britroc_number')
print(somatic_samples_with_no_mutations)

# reformat
somatic_samples_with_no_mutations = somatic_samples_with_no_mutations %>%
  dplyr::select(patient_id) %>% 
  unique() %>% 
  dplyr::mutate(variant_type='dummy', gene_symbol='BRCA1') %>% 
  dplyr::select(patient_id,variant_type,gene_symbol)

somatic_samples_with_mutations = relapse_variants %>% dplyr::select(patient_id,SYMBOL,Consequence)
somatic_samples_with_no_mutations =  somatic_samples_with_no_mutations %>% dplyr::ungroup()
colnames(somatic_samples_with_mutations) = c('patient_id','gene_symbol','variant_type')

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
somatic_variants = somatic_variants %>% tibble::as_tibble() %>% 
  dplyr::group_by(patient_id,gene_symbol) %>%
  dplyr::filter(dplyr::row_number()==1) %>% dplyr::ungroup()

print(somatic_variants %>% dplyr::filter(gene_symbol=='BRCA1'))

# convert the data frame to a matrix
somatic_variants = somatic_variants %>% 
  tidyr::pivot_wider(names_from=patient_id, values_from=variant_type) %>% 
  as.matrix()

row.names(somatic_variants) = somatic_variants[,1]
somatic_variants = somatic_variants[,-1]

# replace NA values
somatic_variants = tidyr::replace_na(somatic_variants, replace = '')

print(somatic_variants)

# map short mutation classes to labels
col = c(
  'frameshift'= 'red',
  'stop_gained'= 'orange',
  'splice_region_SNV'= 'yellow',
  'inframe indel'= 'cyan',
  'missense'= 'blue'
)

alter_fun = list(
  background = ComplexHeatmap::alter_graphic("rect", fill = "#CCCCCC"),   
  frameshift = ComplexHeatmap::alter_graphic("rect", fill = col["frameshift"]),
  stop_gained = ComplexHeatmap::alter_graphic("rect", fill = col["stop_gained"]),
 `inframe indel` = ComplexHeatmap::alter_graphic("rect", fill = col["inframe indel"]),
  splice_region_SNV = ComplexHeatmap::alter_graphic("rect", fill = col["splice_region_SNV"]),
  missense = ComplexHeatmap::alter_graphic("rect", fill = col["missense"])
)

pdf('relapse_britroc_oncoprint.pdf')

ComplexHeatmap::oncoPrint(
  mat=somatic_variants,
  alter_fun = alter_fun,
  row_order = c('TP53','BRCA1','BRCA2','FANCM')
#  heatmap_legend_param = list(labels=c('frameshift','stop gained','inframe indel','splice region SNV','missense'),
#                              title='Alterations')
)

dev.off()

print(relapse_variants)
#print(sequenced_samples, n=Inf)
