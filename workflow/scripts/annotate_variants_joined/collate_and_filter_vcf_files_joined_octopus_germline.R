# A simple script within a larger snakemake workflow of collating and filtering variant calls outputted by a variant calling algorithm

library(DBI)
library(RPostgres)
library(magrittr)

# ensure that the script reads from the users .Renviron text file
readRenviron('~/.Renviron')

vep_files = snakemake@input$vcf_files %>% unlist

print(vep_files %>% length)

print(vep_files)

patient_names = stringr::str_extract(string=vep_files, pattern='[0-9]+.vcf$') %>%
	stringr::str_extract('[0-9]+')

print(patient_names)

# a rudimentary helper function to add the sample IDs to the annotation output table
mutate_x_y = function(x,y) {
  x = x %>% 
	dplyr::rename(
		`sample_1`=dplyr::matches('JBLAB-[0-9]+$'),
		`sample_2`=dplyr::matches('JBLAB-[0-9]+_d$')
		)
  
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
  col_names=TRUE
    ) %>%
  purrr::map2_dfr(.y=patient_names, .f=mutate_x_y ) %>%
  dplyr::select(patient_id, everything())

annotations = annotations %>% unique()

get_vep_variant_format = function (row_index) {

	x = annotations[row_index,]

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

	return(x)
}

vcf = purrr::map_dfr(.x=1:dim(annotations)[1], .f=get_vep_variant_format)

print(vcf)

#collated_vep_file = readr::read_tsv(snakemake@input[['collated_vep_file']])
#print(collated_vep_file)
#collated_vep_file$patient_id = as.character(collated_vep_file$patient_id)


#vcf = 
#	dplyr::semi_join(
#		vcf,
#		collated_vep_file,
#		by=c('patient_id', 'vep_format'='#Uploaded_variation')
#	)
#

print(vcf)




readr::write_tsv(vcf, snakemake@output[[1]])

# editing for MTBP

vcf_MTBP = vcf %>% dplyr::select(`#CHROM`,POS,ID,REF,ALT,QUAL,FILTER,INFO)
vcf_MTBP = vcf_MTBP %>% dplyr::mutate(QUAL=5000, INFO='foo')
vcf_MTBP = vcf_MTBP %>% unique()

# change how indels are indexed to match what is done in the ampliconseq pipeline
vcf_MTBP$POS =
        dplyr::if_else(
                nchar(vcf_MTBP$REF) > 1 | nchar(vcf_MTBP$ALT) > 1,
                vcf_MTBP$POS - 1,
                vcf_MTBP$POS
        )

readr::write_tsv(vcf_MTBP, snakemake@output[[2]], col_names=TRUE)

quit()



print(annotations %>% dplyr::filter(nchar(Allele) != 1) %>% dplyr::select('#Uploaded_variation',Location,Allele), n=100)

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
		      'HGVSc','HGVSp','HGVSg','HGVS_OFFSET'
			)

new_columns = purrr::map_dfc(.x=new_column_names, .f=extract_column)
colnames(new_columns) = new_column_names

annotations = cbind(annotations, new_columns) %>% as_tibble()
annotations$HGVS_united = dplyr::coalesce(annotations$HGVSg, annotations$HGVSc, annotations$HGVSp)

annotations %>% dplyr::filter(is.na(HGVS_united)) %>% print()
quit()

#annotations = annotations %>% dplyr::filter(!grepl(':n.',HGVS_united)) # these are non-coding transcripts for non-target genes for our panels
#annotations = annotations %>% dplyr::filter(!grepl('dup',HGVS_united)) # this looks like it is a bug with MTBP and will need to be fixed at some point
#annotations = annotations %>% dplyr::filter(!grepl('c.\\*',HGVS_united)) # removes three prime utr variants which MTBP is unable to process

annotations = annotations %>% dplyr::filter(
  CANONICAL == 'YES'
)

#annotations = annotations %>% dplyr::filter(
#  !Consequence %in% c(
#    'synonymous_variant'
#  )
#)


#annotations = annotations %>% dplyr::filter(
#  !grepl('benign',CLIN_SIG)
#)

annotations$HGVSp = annotations$HGVSp %>% stringr::str_remove('ENSP[0-9]+\\.[0-9]:')

#annotations %>% dplyr::filter(grepl('+6T',HGVSc)) %>% .$HGVSc %>% print()
#annotations %>% dplyr::filter(!is.na(HGVS_OFFSET)) %>% .$HGVSp %>% print()

readr::write_tsv(annotations %>% dplyr::select(HGVS_united) %>% unique(), col_names=FALSE, file='tmp.txt')
quit()

annotations %>% .$HGVSp %>% unique()

annotations = annotations %>%
	dplyr::filter(
		HGVSp %in% c(
			'p.Lys654SerfsTer47',
			'p.Asn1002ThrfsTer22',
			'p.Lys894ThrfsTer8',
			'p.Glu881Ter',
			'p.Asp825GlufsTer21',
			'p.Glu732Ter',
			'p.His692MetfsTer9',
			'p.Val299ArgfsTer4',
			'p.Lys583AsnfsTer3',
			'p.Arg1203Ter',
			'p.Glu23ValfsTer17',
			'p.Gly1077AlafsTer8',
			'p.Ser1253ArgfsTer10',
			'p.Ala1708Glu',
			'p.Ala1729Glu',
			'p.Ile1465Ter',
			'p.Ile1486Ter',
			'p.Val1736Ala',
			'p.Val1757Ala',
			'p.Tyr1625Ter',
			'p.Tyr1646Ter',
			'p.Ala1623Gly',
			'p.Ala1644Gly',
			'p.Trp1508Ter',
			'p.Trp1529Ter',
			'p.Ser1503Ter',
			'p.Ser1524Ter',
			'p.Pro9GlnfsTer16',
			'p.Glu1035Ter',
			'p.Glu1493ValfsTer10',
			'p.Glu1571GlyfsTer3',
			'p.Asn1784HisfsTer2',
			'p.Ser1982ArgfsTer22',
			'p.Leu2092ProfsTer7',
			'p.Gln2829Ter',
			'p.Thr3033LeufsTer29',
			'p.Ile605TyrfsTer9',
			'p.Asp1469LysfsTer11',
			'p.Gln397LeufsTer25',
			'p.Val220IlefsTer4',
			'p.Ile504SerfsTer22',
			'p.Asn1039IlefsTer2',
			'p.Arg47Ter',
			'p.Gln133Ter',
			'p.Trp268Ter',
			'p.Trp288Ter',
			'p.Lys254ArgfsTer19',
			'p.Arg264Trp'
		) | HGVSc %in% c('ENST00000471181.2:c.302-2del','ENST00000471181.2:c.4548-2A>G','ENST00000471181.2:c.5049+6T>G')
	)

annotations %>% .$HGVSp %>% unique()
annotations %>% .$HGVS_united %>% unique()
#quit()

readr::write_tsv(annotations, path=snakemake@output[[1]], append=FALSE) #'tmp_annotations_joined_archival.tsv'

quit()
#######






annotations %>% dim() %>% print()

## begin to filter mutations for predicted pathogenicity

# remove TP53 and BRCA1 artifacts
annotations = annotations %>% dplyr::filter(!`#Uploaded_variation` %in% c('chr17_7573010_T/G','chr17_7574036_G/A'))
annotations = annotations %>% dplyr::filter(!`#Uploaded_variation` %in% c('chr17_41231352_T/C'))

# filter columns
#annotations = annotations %>% dplyr::filter(
 # !BIOTYPE %in% c('non_stop_decay','nonsense_mediated_decay','processed_pseudogene')
#)

annotations$gnomAD_AF = as.numeric(annotations$gnomAD_AF)
#annotations$AF = as.numeric(annotations$AF)

#annotations = annotations %>% dplyr::filter(
#  gnomAD_AF < 0.001 | is.na(gnomAD_AF)
#)

# 1k genome project
#annotations = annotations %>% dplyr::filter(
#   AF < 0.001  | is.na(AF)
#)

#annotations = annotations %>% dplyr::filter(
# ( !(grepl('tolerated',SIFT) & is.na(CLIN_SIG)) )
#)

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

#annotations = annotations %>% dplyr::filter(
#  !Consequence %in% c(
#    'downstream_gene_variant','5_prime_UTR_variant','3_prime_UTR_variant',
#    'upstream_gene_variant', 'non_coding_transcript_exon_variant',
#    'intron_variant,non_coding_transcript_variant', 'intron_variant',
#    'synonymous_variant'
#  )
#)


# remove transcript isoform duplicates

# this doesn't make a change if we already filter for 'CANONICAL=TRUE'

#annotations = annotations %>%
#  dplyr::group_by(sample_id,Location,Allele) %>%
#  dplyr::filter(dplyr::row_number()==1) %>% dplyr::ungroup()

print(annotations %>% dplyr::select(patient_id,`#Uploaded_variation`,Allele), n=Inf)

vcf = dplyr::inner_join(vcf,annotations, by=c('vep_format'='#Uploaded_variation', 'patient_id'))

readr::write_tsv(annotations, path=snakemake@output[[1]], append=FALSE) #'tmp_annotations_joined_archival.tsv'

annotations %>% group_by(SYMBOL,`#Uploaded_variation`) %>% summarise(n=n()) %>% arrange(n) %>% print(n=Inf)
annotations %>% group_by(patient_id,SYMBOL) %>% summarise(n=n()) %>% dplyr::ungroup() %>% dplyr::group_by(SYMBOL) %>% 
	summarise(n=n()) %>% dplyr::arrange(n) %>% print(n=Inf)

quit()
