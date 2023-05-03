# A simple script within a larger snakemake workflow of collating and filtering variant calls outputted by a variant calling algorithm

library(magrittr)

# ensure that the script reads from the users .Renviron text file
readRenviron('~/.Renviron')

vep_files = snakemake@input$vep_files %>% unlist

patient_names = stringr::str_extract(string=vep_files, pattern='[0-9]+(.filtered)?.vep.vcf$') %>%
	stringr::str_extract('[0-9]+')

# a rudimentary helper function to add the patient IDs to the annotation output table
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

# filter to ensure tech reps had matching genotypes
vcf = readr::read_tsv(snakemake@input$vcf_file)
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
	} else if (nchar(x$REF) > 1) { #deletion
		x$POS = x$POS + 1
		x$REF_tmp = x$REF %>% stringr::str_remove('[ATGC]')
		x$ALT_tmp = x$ALT %>% stringr::str_replace('[ATGC]', replacement='-')
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

vcf = purrr::map_dfr(.x=1:dim(vcf)[1], .f=get_vep_variant_format)
vcf$patient_id = vcf$patient_id %>% as.character()

annotations = dplyr::inner_join(annotations,vcf, by=c('#Uploaded_variation'='vep_format', 'patient_id'))

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

## begin to filter mutations for predicted pathogenicity

# remove TP53 and BRCA1 artifacts
#annotations = annotations %>% dplyr::filter(!`#Uploaded_variation` %in% c('chr17_7573010_T/G','chr17_7574036_G/A'))
#annotations = annotations %>% dplyr::filter(!`#Uploaded_variation` %in% c('chr17_41231352_T/C'))

# filter columns
annotations$gnomAD_AF = as.numeric(annotations$gnomAD_AF)

annotations = annotations %>% dplyr::filter(
  CANONICAL == 'YES'
)

# this called a splice mutation as low impact so probably not a good filter

annotations = annotations %>% dplyr::filter(
  !Consequence %in% c(
    'downstream_gene_variant', 'upstream_gene_variant'
  )
)

annotations = annotations %>% dplyr::filter(SYMBOL!='RAD51L3-RFFL')

vcf = dplyr::inner_join(vcf,annotations, by=c('vep_format'='#Uploaded_variation', 'patient_id'))

annotations$patient_id = annotations$patient_id %>% as.integer()
annotations = annotations %>% dplyr::arrange(patient_id)
vep_reduced = annotations %>% dplyr::select(patient_id, SYMBOL, CHROM, POS, REF, ALT, Consequence, HGVSp, HGVSc)
vep_reduced$HGVSc = vep_reduced$HGVSc %>% stringr::str_remove('ENST[0-9]+\\.[0-9]+:')
vep_reduced$HGVSp = vep_reduced$HGVSp %>% stringr::str_remove('ENSP[0-9]+\\.[0-9]+:')

vep_reduced = vep_reduced %>% dplyr::select(patient_id, SYMBOL, CHROM, POS, REF, ALT, Consequence, HGVSc, HGVSp) %>% unique()

print(vep_reduced)
print(snakemake@output[['vep_reduced']])

readr::write_tsv(annotations, file=snakemake@output[['vep_output']], append=FALSE) #'tmp_annotations_joined_archival.tsv'
readr::write_tsv(vep_reduced, file=snakemake@output[['vep_reduced']], append=FALSE) #'tmp_annotations_joined_archival.tsv'
readr::write_tsv(vcf, file=snakemake@output[['vcf_output']], append=FALSE) #'tmp_annotations_joined_archival.tsv'

quit()
