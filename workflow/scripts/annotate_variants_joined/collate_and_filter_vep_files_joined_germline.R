# A simple script within a larger snakemake workflow of collating and filtering variant calls outputted by a variant calling algorithm

library(magrittr)
library(DBI)
library(RPostgres)

# ensure that the script reads from the users .Renviron text file
readRenviron('~/.Renviron')

vep_files = snakemake@input$vep_files %>% unlist

print(vep_files %>% length)

print(vep_files)

patient_names = stringr::str_extract(string=vep_files, pattern='[0-9]+.filtered.vep.vcf$') %>%
	stringr::str_extract('[0-9]+')

print(patient_names)

# a rudimentary helper function to add the sample IDs to the annotation output table
mutate_x_y = function(x,y) {
  return(dplyr::mutate(x, fk_britroc_number=y))
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
  dplyr::select(fk_britroc_number, everything())

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
	x$vep_format = paste(x$`#CHROM`, x$POS, x$Alleles, sep='_')
	x$Alleles = NULL

	print(x)

	return(x)
}

print(vcf)

print('foo')

vcf = purrr::map_dfr(.x=1:dim(vcf)[1], .f=get_vep_variant_format)
vcf$patient_id = vcf$patient_id %>% as.character() # snakemake@wildcards[['patient_id']]

print('goo')

print(vcf)

print(annotations)
print(vcf, width=Inf)

#quit()

annotations = dplyr::inner_join(annotations,vcf, by=c('#Uploaded_variation'='vep_format', 'fk_britroc_number'='patient_id'))

print('shoe')

print(annotations)

#quit()

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

annotations = cbind(annotations, new_columns) %>% tibble::as_tibble()

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

print(annotations %>% dplyr::select(fk_britroc_number,`#Uploaded_variation`,Allele), n=Inf)

vcf = dplyr::inner_join(vcf,annotations, by=c('vep_format'='#Uploaded_variation', 'patient_id'='fk_britroc_number'))

readr::write_tsv(annotations, path=snakemake@output[[1]], append=FALSE) #'tmp_annotations_joined_archival.tsv'
#readr::write_tsv(vcf, path='tmp_annotations_joined_archival2.tsv', append=FALSE)

#annotations %>% group_by(SYMBOL,`#Uploaded_variation`) %>% summarise(n=n()) %>% dplyr::arrange(n) %>% print(n=Inf)
#annotations %>% group_by(patient_id,SYMBOL) %>% summarise(n=n()) %>% dplyr::ungroup() %>% dplyr::group_by(SYMBOL) %>% 
#	summarise(n=n()) %>% dplyr::arrange(n) %>% print(n=Inf)

quit()
