library(magrittr)
library(DBI)
library(RPostgres)
library(ComplexHeatmap)

source('~/.Renviron')

# read in tumour samples variants
non_tp53_variants = readr::read_tsv(snakemake@input[['filtered_non_TP53_variants']])

non_tp53_variants_original = non_tp53_variants

# tumour type as classified by the study design
tumour_type = snakemake@wildcards$tumour_sample_type

# Set of annotated TP53 variants 
tp53_variants = readr::read_tsv(snakemake@input[['filtered_TP53_variants_with_MAFs']]) %>% 
		dplyr::filter(type==tumour_type) %>% 
		dplyr::mutate(SYMBOL='TP53')

print(tp53_variants)
print(non_tp53_variants)

# removed variants on the basis of annotations in the MTBP pipeline
if (tumour_type=='archival') {
non_tp53_variants = non_tp53_variants %>% 
	dplyr::filter(`#Uploaded_variation`!='chr17_41231352_T/C') %>% # suspected artifact
	dplyr::filter(!(SYMBOL == 'FANCM' & patient_id==69))
} else {
non_tp53_variants = non_tp53_variants %>% 
	dplyr::filter(`#Uploaded_variation`!='chr17_41231352_T/C') %>% # suspected artifact
	dplyr::filter(!(SYMBOL == 'FANCM' & patient_id==123)) %>%
	dplyr::filter(!(SYMBOL == 'BRCA2' & patient_id==77))
}

# read in germline variant information
germline_variants = readr::read_tsv(snakemake@input[['filtered_germline_variants']]) %>% 
	dplyr::rename(patient_id='fk_britroc_number', Consequence='variant_type', SYMBOL='gene_symbol')
germline_variants = germline_variants %>% dplyr::select(patient_id,Consequence,SYMBOL) %>% dplyr::mutate(type='germline')

# recode germline variant information
germline_variants$Consequence = dplyr::recode(
	germline_variants$Consequence,
	'splice_acceptor'='splice_acceptor_variant',
	'frameshift'='frameshift_variant',
	'nonsynonymous'='missense_variant',
	'stop_gained,splice_region'='stop_gained,splice_region_variant',
	'splice_region,intron'='splice_region_variant,intron_variant',
	'framehift,splice_region'='frameshift_variant'
)

# read in TP53 variant clonality predictions
tp53_variant_clonality = readr::read_tsv(snakemake@input[['clonality_status_of_TP53_variants']])

# join the two TP53 variant tables
tp53_variants_clonal_mutations = dplyr::inner_join(
	tp53_variants, 
	tp53_variant_clonality %>% dplyr::filter(classification=='clonal' | is.na(classification)), 
	by=c('fk_britroc_number','#Uploaded_variation')
) %>% dplyr::select(fk_britroc_number,Consequence,SYMBOL)

# retain non-clonal variants for patients in which a clonal TP53 hasn't been found
tp53_variants_no_clonal_mutation = tp53_variants %>% 
	dplyr::filter(!fk_britroc_number %in% tp53_variants_clonal_mutations$fk_britroc_number) %>%
	dplyr::select(fk_britroc_number,Consequence,SYMBOL)

# bind TP5 variants together
tp53_variants = rbind(tp53_variants_clonal_mutations, tp53_variants_no_clonal_mutation)

# reformat variant information
non_tp53_variants = non_tp53_variants %>% dplyr::select(patient_id,Consequence,SYMBOL) %>% dplyr::mutate(type='somatic')
tp53_variants = tp53_variants %>% dplyr::rename(patient_id=fk_britroc_number) %>% dplyr::mutate(type='somatic')

# filter variants by type and by gene symbol

if (tumour_type=='archival') {
	non_tp53_variants = non_tp53_variants %>% dplyr::filter(Consequence %in% c('frameshift_variant','stop_gained'))
	non_tp53_variants = non_tp53_variants %>% dplyr::filter(SYMBOL %in% c('BRCA1','BRCA2','FANCM','BARD1','RAD51B','RAD51C','RAD51D','PALB2'))
} else  {
	non_tp53_variants = non_tp53_variants %>% dplyr::filter(Consequence %in% c('frameshift_variant','stop_gained'))
	non_tp53_variants = non_tp53_variants %>% dplyr::filter(SYMBOL %in% c('BRCA1','BRCA2','FANCM','BARD1','RAD51B','RAD51C','RAD51D','PALB2'))
}

# bind all variants together
all_variants = rbind(tp53_variants, non_tp53_variants, germline_variants)

# factorise variant consequence column
all_variants$Consequence = factor(all_variants$Consequence,
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

# reformat column
all_variants = all_variants %>% dplyr::arrange(patient_id,Consequence)

# TODO: formalise this
# add variants that were previously filtered due to an overly harsh MAF threshold

if (tumour_type=='archival') {
	all_variants = all_variants %>% tibble::add_row(SYMBOL='BRCA2', patient_id=39, Consequence='frameshift', type='somatic')
	all_variants = all_variants %>% tibble::add_row(SYMBOL='BRCA2', patient_id=141, Consequence='frameshift', type='somatic')
} else if (tumour_type=='relapse') {
	all_variants = all_variants %>% tibble::add_row(SYMBOL='FANCM', patient_id=176, Consequence='frameshift', type='somatic')
}

# retrieve sequencing metadata
britroc_con <- dbConnect(RPostgres::Postgres(),
                         dbname='britroc1',
                         host='jblab-db.cri.camres.org',
                         port = 5432,
                         user = jblab_db_username,
                         password = jblab_db_password
)
clarity_con <- dbConnect(RPostgres::Postgres(),
                         dbname='clarity',
                         host='jblab-db.cri.camres.org',
                         port = 5432,
                         user = jblab_db_username,
                         password = jblab_db_password
)

# read in DNA sample information
samples = dbReadTable(britroc_con, 'sample')

slx_library = dbReadTable(britroc_con, 'slx_library') %>%
       dplyr::filter(fk_slx!='SLX-13716') %>% # no data for this SLX
       dplyr::filter(grepl('AA',fk_experiment)) # select for tam-seq experiments only

slx_clarity = dbReadTable(clarity_con, 'slx')

slx_library = dplyr::semi_join(slx_library, slx_clarity, by=c('fk_slx'='name')) # ensure slx is actually in clarity

experiments = dbReadTable(britroc_con, 'experiment') %>% dplyr::filter(fk_amplicon_panel %in% c(6,28)) # only allow panel 6 or panel 28 amplicon panels
slx_library = slx_library %>% dplyr::semi_join(experiments, by=c('fk_experiment'='name'))

# retrieve analysis type
sequenced_samples = dplyr::semi_join(samples, slx_library, by=c('name'='fk_sample'))

# remove non-HGSOC samples
non_hgsoc_samples = c('JBLAB-4114','JBLAB-4916','IM_249','IM_250','IM_234','IM_235','IM_236','IM_237','JBLAB-4271','IM_420',
			'IM_262','JBLAB-4922','JBLAB-4923','IM_303','IM_290','IM_43','IM_293','IM_307','IM_308','IM_309','IM_424',
			'IM_302','IM_303','IM_304','IM_305',
			'JBLAB-19320','IM_61','IM_62','IM_63','IM_397','IM_302','IM_98','JBLAB-4210','IM_147','JBLAB-4216','IM_44')

# remove samples with poor sequencing QC
samples_with_no_good_sequencing = c('IM_144','IM_435','IM_436','IM_158','IM_296','IM_373','IM_154','IM_297','IM_365','IM_432','IM_429','IM_368','IM_441') 

# remove samples with very low purity
samples_with_very_low_purity = c('IM_1','IM_2','IM_3','IM_4','IM_20','IM_26','IM_27','IM_69','IM_86','IM_90','IM_93','IM_94','IM_173','IM_177','IM_179','IM_200',
				'IM_241','IM_242','IM_417','IM_418','IM_419','IM_420','IM_221','IM_264','IM_329','IM_289','IM_308','IM_309',
				'IM_338','IM_339','IM_340','IM_341','IM_342','IM_432','IM_372','IM_272','IM_392')

samples_which_may_have_been_mislabelled = c('JBLAB-4119','IM_92','IM_95','IM_116','IM_413')

# bind the poor samples together
bad_samples = c(non_hgsoc_samples, samples_with_no_good_sequencing, samples_with_very_low_purity, samples_which_may_have_been_mislabelled)

# remove samples from downstream analyses
sequenced_samples = sequenced_samples %>% dplyr::filter(!name %in% bad_samples)

# group sequenced samples by patient, and only retain patients with at least one germline, relapse and archival sample
sequenced_samples = sequenced_samples %>%
       dplyr::group_by(fk_britroc_number,type) %>%
       dplyr::summarise(n=dplyr::n()) %>%
       tidyr::pivot_wider(names_from=type, values_from=n) %>%
       dplyr::filter(!is.na(archival) & !is.na(relapse) & !is.na(germline)) # https://stackoverflow.com/questions/27197617/filter-data-frame-by-character-column-name-in-dplyr

# only retain variants in patients with at least one germline, relapse and archival sample
all_variants = all_variants %>% dplyr::filter(patient_id %in% sequenced_samples$fk_britroc_number)

# reformat
all_variants =
  all_variants %>% dplyr::select(patient_id,SYMBOL,Consequence,type) %>%
  unique()

# remove synonymous and intron variants
all_variants = 
  all_variants %>% dplyr::filter(!Consequence %in% c('intron_variant','synonymous_variant')) %>%
  unique()

# this information is needed so that the oncoprint shows information even for patients without a mutation in any one gene
somatic_samples_with_no_mutations = dplyr::anti_join(
  sequenced_samples, all_variants, by=c('fk_britroc_number'='patient_id')
)

somatic_samples_with_no_mutations = somatic_samples_with_no_mutations %>% dplyr::rename(patient_id='fk_britroc_number')

# reformat
somatic_samples_with_no_mutations = somatic_samples_with_no_mutations %>%
  dplyr::select(patient_id) %>% 
  unique() %>% 
  dplyr::mutate(variant_type='dummy', gene_symbol='BRCA1') %>% 
  dplyr::select(patient_id,variant_type,gene_symbol)

somatic_samples_with_mutations = all_variants %>% dplyr::select(patient_id,SYMBOL,Consequence,type)
somatic_samples_with_no_mutations =  somatic_samples_with_no_mutations %>% dplyr::ungroup() %>% dplyr::mutate(mutation_type='somatic')
colnames(somatic_samples_with_mutations) = c('patient_id','gene_symbol','variant_type','mutation_type')

print(somatic_samples_with_mutations)
print(somatic_samples_with_no_mutations)

print('foo')
somatic_variants = rbind(somatic_samples_with_mutations, somatic_samples_with_no_mutations)
print('goo')

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
  dplyr::filter(dplyr::row_number()==1) %>% 
  dplyr::ungroup() %>%
  dplyr::arrange(patient_id)

print(somatic_variants)

#non_tp53_variants = non_tp53_variants_original %>% dplyr::inner_join(somatic_variants, by=c('patient_id', 'SYMBOL'='gene_symbol'))

#print(non_tp53_variants)

readr::write_tsv(somatic_variants %>% dplyr::filter(mutation_type=='somatic'), 'BriTROC-1_relapse_somatic_mutations.tsv', append=FALSE)
#readr::write_tsv(somatic_variants %>% dplyr::filter(mutation_type=='germline'), 'BriTROC-1_germline_mutations.tsv', append=FALSE)

quit()

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

print(somatic_variants)
dim(somatic_variants)

dummy = matrix('', nrow=6, ncol=108)
row.names(dummy) = c('BARD1','RAD51B','RAD51C','RAD51D','BRIP1','PALB2')

dim(somatic_variants) %>% print()
dim(dummy) %>% print()

x = rbind(somatic_variants, dummy)

z = x[1,]
y = x[1,] %>% tibble::as_tibble() %>% dplyr::mutate(patient_id=names(z))

print(y %>% dplyr::filter(value==''))

somatic_variants=x

# generate a platinum sensitivity partition
heatmap_patients = colnames(somatic_variants)
heatmap_table = tibble::tibble(fk_britroc_number=as.integer(heatmap_patients))
print(heatmap_table)

patients = dbReadTable(britroc_con, 'patients')
print(patients)

heatmap_table = dplyr::inner_join(
  heatmap_table, patients,
  by=c('fk_britroc_number'='britroc_number')
) %>% dplyr::select('fk_britroc_number','pt_sensitivity_at_reg')

heatmap_table$pt_sensitivity_at_reg =
  dplyr::recode(
    heatmap_table$pt_sensitivity_at_reg,
    'sensitive'='platinum sensitive',
    'resistant'='platinum resistant'
  )

print(heatmap_table)

pdf(snakemake@output[['germline_and_somatic_oncoprint']])

if (tumour_type=='archival') {
	genes_with_variants = c('TP53','BRCA1','BRCA2','FANCM','BARD1','RAD51B','RAD51C','RAD51D','BRIP1','PALB2')
} else if (tumour_type=='relapse') {
	genes_with_variants = c('TP53','BRCA1','BRCA2','FANCM','BARD1','RAD51B','RAD51C','RAD51D','BRIP1','PALB2')
}

#ComplexHeatmap::oncoPrint(
#  mat=somatic_variants,
#  alter_fun = alter_fun,
#  row_order = genes_with_variants
#)

#ha = HeatmapAnnotation(show_legend = c(FALSE))

ComplexHeatmap::oncoPrint(
  mat=somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum resistant'],
  alter_fun = alter_fun,
  show_heatmap_legend=FALSE,
  row_labels=c('','','','','','','','','',''),
  column_title='platinum resistant',
  row_order = genes_with_variants
) +
ComplexHeatmap::oncoPrint(
  mat=somatic_variants[,heatmap_table$pt_sensitivity_at_reg=='platinum sensitive'],
  alter_fun = alter_fun,
  column_title='platinum sensitive',
 # heatmap_legend_param = list(labels=c('frameshift','stop gained','splice region SNV','missense'),
#                              title='Alterations'),
  row_order = genes_with_variants
)

dev.off()
