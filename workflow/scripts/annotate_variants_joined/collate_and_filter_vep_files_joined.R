# A simple script within a larger snakemake workflow of collating and filtering variant calls outputted by a variant calling algorithm

library(tidyverse)
library(DBI)
library(RPostgres)

# ensure that the script reads from the users .Renviron text file
readRenviron('~/.Renviron')

vep_files = snakemake@input$vep_files %>% unlist

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
      `#Uploaded_variation` = 'c',
      Location = 'c',
      cDNA_position='c',
      CDS_position='c',
      Protein_position='c',
      Allele='c',
      Amino_acids='c'
    )
)  %>%
  purrr::map2_dfr(.y=patient_names, .f=mutate_x_y ) %>%
  dplyr::arrange(Gene) %>%
  dplyr::select(patient_id, everything())

annotations = annotations %>% unique()
print(annotations %>% dplyr::filter(nchar(Allele) != 1) %>% dplyr::select('#Uploaded_variation',Location,Allele), n=100)


# filter to ensure tech reps had matching genotypes
vcf = readr::read_tsv(snakemake@input$vcf_file)
#vcf = vcf %>% tidyr::unite('#CHROM',POS, col='Location', sep=':')
get_vep_variant_format = function (row_index) {

	x = vcf[row_index,]

	if (grepl('\\*',x$ALT) & nchar(x$REF) == 1) {
		print('complex multi-allelic variant')
		x$ALT_tmp = x$ALT %>% stringr::str_replace(',', replacement='/')
		x$REF_tmp = x$REF
	}
	else if(grepl('\\*',x$ALT)) {
		print('complex multi-allelic variant')
		x$POS = x$POS + 1
		x$REF_tmp = x$REF %>% stringr::str_remove('[ATGC]')
		x$ALT_tmp = x$ALT %>% stringr::str_replace('[ATGC]', replacement='-')
		x$ALT_tmp = x$ALT_tmp %>% stringr::str_replace(',', replacement='/')	
	}
	else if (nchar(x$REF) > 1 & nchar(x$ALT) > 1) { 
		x$REF_tmp = x$REF
		x$ALT_tmp = x$ALT
	}
	else if (nchar(x$ALT) > 1) { #insertion
		x$POS = x$POS + 1
		x$REF_tmp = x$REF %>% stringr::str_replace('[ATGC]', replacement='-')
		x$ALT_tmp = x$ALT %>% stringr::str_remove('[ATGC]')
		#print(x$REF)
		#print(x$REF_tmp)
		#print(x$ALT)
		#print(x$ALT_tmp)
	} else if (nchar(x$REF) > 1) { #deletion
		x$POS = x$POS + 1
		x$REF_tmp = x$REF %>% stringr::str_remove('[ATGC]')
		x$ALT_tmp = x$ALT %>% stringr::str_replace('[ATGC]', replacement='-')
		#print(x$REF)
		#print(x$REF_tmp)
		#print(x$ALT)
		#print(x$ALT_tmp)
	}
	else {
		x$REF_tmp = x$REF
		x$ALT_tmp = x$ALT
	}

	x = x %>% tidyr::unite(REF_tmp, ALT_tmp, col='Alleles', sep='/', remove=TRUE)
	x$vep_format = paste(x$CHROM, x$POS, x$Alleles, sep='_')
	x$Alleles = NULL

	return(x)
}

print(vcf)

print('foo')

vcf = purrr::map_dfr(.x=1:dim(vcf)[1], .f=get_vep_variant_format)
vcf$patient_id = vcf$patient_id %>% as.character()

print('goo')

print(vcf)

print(annotations)

annotations = dplyr::inner_join(annotations,vcf, by=c('#Uploaded_variation'='vep_format', 'patient_id'))

print('shoe')

print(annotations)

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
                     'ENSP','UNIPARC','EXON','INTRON','SWISSPROT','TREMBL','DOMAINS','FLAGS','gnomAD_AF','AF',
		      'HGVSc','HGVSp'
			)

new_columns = purrr::map_dfc(.x=new_column_names, .f=extract_column)
colnames(new_columns) = new_column_names

annotations = cbind(annotations, new_columns) %>% as_tibble()

annotations %>% dim() %>% print()

## begin to filter mutations for predicted pathogenicity

# remove TP53 and BRCA1 artifacts
annotations = annotations %>% dplyr::filter(!`#Uploaded_variation` %in% c('chr17_7573010_T/G','chr17_7574036_G/A'))
annotations = annotations %>% dplyr::filter(!`#Uploaded_variation` %in% c('chr17_41231352_T/C'))

# filter columns
annotations = annotations %>% dplyr::filter(
  !BIOTYPE %in% c('non_stop_decay','nonsense_mediated_decay','processed_pseudogene')
)

annotations$gnomAD_AF = as.numeric(annotations$gnomAD_AF)
#annotations$AF = as.numeric(annotations$AF)

annotations = annotations %>% dplyr::filter(
  gnomAD_AF < 0.001 | is.na(gnomAD_AF)
)

# 1k genome project
#annotations = annotations %>% dplyr::filter(
#   AF < 0.001  | is.na(AF)
#)

annotations = annotations %>% dplyr::filter(
  !grepl('benign',CLIN_SIG)
)

annotations = annotations %>% dplyr::filter(
 ( !(grepl('tolerated',SIFT) & is.na(CLIN_SIG)) )
)

annotations = annotations %>% dplyr::filter(
  ( !(grepl('benign',PolyPhen) & is.na(CLIN_SIG)) )
)

#annotations = annotations %>% dplyr::filter(
#  !( grepl('possibly',PolyPhen) & grepl('low_confidence',SIFT))
#)

#grepl('possibly',annotations$PolyPhen)
#grepl('low_confidence',annotations$SIFT) 
#c(grepl('possibly',annotations$PolyPhen)) & c(grepl('low_confidence',annotations$SIFT))

#annotations = annotations %>% dplyr::filter(
#  SYMBOL %in% c('BARD1','BRCA1','BRCA2','BRIP1','FANCM','PALB2','RAD51B','RAD51C','RAD51D')
#)

annotations = annotations %>% dplyr::filter(
  CANONICAL == 'YES'
)

# this called a splice mutation as low impact so probably not a good filter

#annotations = annotations %>% dplyr::filter(
#  IMPACT != 'LOW'
#)

#annotations = annotations %>% dplyr::filter(
#  !Consequence %in% c(
#    '3_prime_UTR_variant','5_prime_UTR_variant','downstream_gene_variant','intron_variant',
#    'upstream_gene_variant', 'non_coding_transcript_exon_variant',
#    'intron_variant,non_coding_transcript_variant'
#  )
#)

annotations = annotations %>% dplyr::filter(
  !Consequence %in% c(
    'downstream_gene_variant','5_prime_UTR_variant','3_prime_UTR_variant',
    'upstream_gene_variant', 'non_coding_transcript_exon_variant',
    'intron_variant,non_coding_transcript_variant', 'intron_variant',
    'synonymous_variant'
  )
)


# remove transcript isoform duplicates

# this doesn't make a change if we already filter for 'CANONICAL=TRUE'

#annotations = annotations %>%
#  dplyr::group_by(sample_id,Location,Allele) %>%
#  dplyr::filter(dplyr::row_number()==1) %>% dplyr::ungroup()

print(annotations %>% dplyr::select(patient_id,`#Uploaded_variation`,Allele), n=Inf)

vcf = dplyr::inner_join(vcf,annotations, by=c('vep_format'='#Uploaded_variation', 'patient_id'))

readr::write_tsv(annotations, path=snakemake@output[[1]], append=FALSE) #'tmp_annotations_joined_archival.tsv'
readr::write_tsv(vcf, path='tmp_annotations_joined_archival2.tsv', append=FALSE)

annotations %>% group_by(SYMBOL,`#Uploaded_variation`) %>% summarise(n=n()) %>% arrange(n) %>% print(n=Inf)
annotations %>% group_by(patient_id,SYMBOL) %>% summarise(n=n()) %>% dplyr::ungroup() %>% dplyr::group_by(SYMBOL) %>% 
	summarise(n=n()) %>% dplyr::arrange(n) %>% print(n=Inf)

quit()

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

annotations = dplyr::inner_join(annotations,samples, by=c('sample_id'='name'))

### for this tumour type calculate the number of people with paired samples ####

#print(vep_files %>% length)

vep_sample_list = stringr::str_extract(string=vep_files, pattern='(IM_[0-9]+)|JBLAB-[0-9]+')
non_hgsoc_samples = c('JBLAB-4114','JBLAB-4916','IM_249','IM_250','IM_234','IM_235','IM_236','IM_237','JBLAB-4271','IM_420',
			'IM_262','JBLAB-4922','JBLAB-4923','IM_303','IM_290','IM_43','IM_293','IM_307','IM_308','IM_309','IM_424',
			'IM_302','IM_303','IM_304','IM_305',
			'JBLAB-19320','IM_61','IM_62','IM_63','IM_397','IM_302','IM_98','JBLAB-4210','IM_147','JBLAB-4216','IM_44')


# IM_44 was not sequenced for sWGS but very likely non-HGSOC as IM_43 was identified as non-HGSOC

samples_with_no_good_sequencing = c('IM_144','IM_435','IM_436','IM_158','IM_296','IM_373','IM_154','IM_297','IM_365','IM_432','IM_429','IM_368','IM_441') 

samples_with_very_low_purity = c('IM_1','IM_2','IM_3','IM_4','IM_20','IM_26','IM_27','IM_69','IM_86','IM_90','IM_93','IM_94','IM_173','IM_177','IM_179','IM_200',
				'IM_241','IM_242','IM_417','IM_418','IM_419','IM_420','IM_221','IM_264','IM_329','IM_289','IM_308','IM_309',
				'IM_338','IM_339','IM_340','IM_341','IM_342','IM_432','IM_372','IM_272','IM_392')

vep_sample_list = vep_sample_list[!vep_sample_list %in% non_hgsoc_samples]
vep_sample_list = vep_sample_list[!vep_sample_list %in% samples_with_no_good_sequencing]
vep_sample_list = vep_sample_list[!vep_sample_list %in% samples_with_very_low_purity]

annotations = annotations %>% filter(sample_id %in% vep_sample_list) # evidence that we may be filtering genuine variants
x = annotations %>% group_by(fk_britroc_number,SYMBOL) %>% summarise(n=n()) %>% ungroup %>% group_by(SYMBOL) %>% summarise(n=n()) %>% arrange(n)

num_patients_with_paired_samples = samples %>% dplyr::filter(name %in% vep_sample_list) %>% .$fk_britroc_number %>% as.factor %>% levels %>% length()
patients_with_paired_samples = samples %>% dplyr::filter(name %in% vep_sample_list) %>% .$fk_britroc_number %>% as.factor %>% levels()

#samples %>% dplyr::filter(name %in% vep_sample_list) %>% print()

print('num patients with sequencing')
print(num_patients_with_paired_samples)

#non_hgsoc_patients = c(10,19,27,111,119,120,144,157,164,170,176,183,197,243,258,274)

# determine samples without a variant after the filtering which occurs in this script

patients_sequenced = samples %>% dplyr::filter(name %in% vep_sample_list)
print(patients_with_paired_samples[!patients_with_paired_samples %in% annotations$fk_britroc_number] %>% unique())
#print(vep_sample_list[!vep_sample_list %in% annotations$sample_id])

#num_patients_with_paired_samples = 141

#slx_library = dbReadTable(britroc_con, 'slx_library') %>% 
#	dplyr::filter(fk_slx!='SLX-13716') %>% # no data for this SLX
#	dplyr::filter(grepl('AA',fk_experiment)) # select for tam-seq experiments only

#slx_clarity = dbReadTable(clarity_con, 'slx')

#slx_library = dplyr::semi_join(slx_library, slx_clarity, by=c('fk_slx'='name')) # ensure slx is actually in clarity

#experiments = dbReadTable(britroc_con, 'experiment') %>% dplyr::filter(fk_amplicon_panel %in% c(6,28)) # only allow panel 6 or panel 28 amplicon panels
#slx_library = slx_library %>% dplyr::semi_join(experiments, by=c('fk_experiment'='name'))

# retrieve analysis type
#analysis_type = stringr::str_extract(string=snakemake@output[[1]], pattern='(archival|relapse)')
#print(analysis_type)
	
#sequenced_samples = dplyr::semi_join(samples, slx_library, by=c('name'='fk_sample')) 
#sequenced_samples = sequenced_samples %>% 
#	dplyr::group_by(fk_britroc_number,type) %>% 
#	dplyr::summarise(n=n()) %>% 
#	tidyr::pivot_wider(names_from=type, values_from=n) %>%
#	dplyr::filter(!is.na(!!as.symbol(analysis_type)) & !is.na(germline)) # https://stackoverflow.com/questions/27197617/filter-data-frame-by-character-column-name-in-dplyr

#num_patients_with_paired_samples = sequenced_samples %>% .$fk_britroc_number %>% as.factor %>% levels %>% length()
#print('number of patients with paired samples for this analysis type:')
#print(num_patients_with_paired_samples)

##################################################################################

x = x %>% dplyr::mutate(prop_patients_with_mutation=n/num_patients_with_paired_samples)
print(x)

print('BRCA1 patients')

annotations %>% dplyr::filter(SYMBOL=='BRCA1') %>% .$fk_britroc_number %>% unique()

print('BRCA2 patients')

annotations %>% dplyr::filter(SYMBOL=='BRCA2') %>% .$fk_britroc_number %>% unique()

print('Total BRCA somatic mutation rate - patient level')

z = annotations %>% dplyr::filter(SYMBOL %in% c('BRCA1','BRCA2')) %>% .$fk_britroc_number %>% unique() %>% length
print(z / num_patients_with_paired_samples)

# add records for samples without any mutations
samples_with_no_mutations = vep_sample_list[!vep_sample_list %in% annotations$sample_id]
samples_with_no_mutations_annotations = annotations[1:length(samples_with_no_mutations),]

samples_with_no_mutations_annotations = samples_with_no_mutations_annotations[NA,]
samples_with_no_mutations_annotations$sample_id = samples_with_no_mutations
samples_with_no_mutations_annotations = samples_with_no_mutations_annotations %>% dplyr::select(-type,-fk_britroc_number)
samples_with_no_mutations_annotations = dplyr::inner_join(samples_with_no_mutations_annotations, samples, by=c('sample_id'='name'))

print(samples_with_no_mutations_annotations)

annotations = rbind(annotations,samples_with_no_mutations_annotations)

write_tsv(annotations, snakemake@output[[1]])

