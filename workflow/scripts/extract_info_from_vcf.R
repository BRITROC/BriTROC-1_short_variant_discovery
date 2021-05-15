library(magrittr)

#vcf = readr::read_tsv('results/tumour_sample_vcfs_octopus/IM_71.filtered2.vcf', comment='##')
#vcf = vcf %>% dplyr::select('#CHROM','POS','REF','ALT','QUAL','FORMAT', starts_with('IM'), starts_with('JBLAB'))

sample_names = snakemake@input %>% stringr::str_extract(pattern='(IM_[0-9]+|JBLAB-[0-9]+)')

list_of_vcf_dfs = purrr::map(
	.x=snakemake@input, 
	.f=readr::read_tsv, 
	comment='##',
	guess_max=5000,
	col_types = 
	readr::cols(
		POS = 'i',
		QUAL='d'
		)
	)

rename_vcf_sample_columns = function(vcf_df) {
	vcf_df = vcf_df %>% dplyr::rename_at(dplyr::vars(names(.) %>% tail(2)), dplyr::funs(paste(c('tech_rep_1','tech_rep_2'))))
	return(vcf_df)
}

mutate_x_y = function(x,y) {
  return(dplyr::mutate(x, sample_id=y))
}

vcf = purrr::map(
	.x=list_of_vcf_dfs, 
	.f=rename_vcf_sample_columns
	) %>% 
	purrr::map2_dfr(.y=sample_names, .f=mutate_x_y)

print(vcf)

#print(vcf)

#vcf$variant_type = ifelse(nchar(vcf$REF) == 1 & nchar(vcf$ALT) == 1, 'SNV', NA)
#vcf$variant_type = ifelse(nchar(vcf$REF) > 1, 'Deletion', vcf$variant_type)
#vcf$variant_type = ifelse(nchar(vcf$ALT) > 1, 'Insertion', vcf$variant_type)

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

	return(x)
}

vcf = purrr::map_dfr(.x=1:dim(vcf)[1], .f=get_vep_variant_format)

vcf$FORMAT = vcf$FORMAT %>% stringr::str_split(':')
vcf$tech_rep_1 = vcf$tech_rep_1 %>% stringr::str_split(':')
vcf$tech_rep_2 = vcf$tech_rep_2 %>% stringr::str_split(':')

get_AFs = function(row_index) {

	vcf = vcf[row_index,]	

	vcf$FORMAT[[1]] = vcf$FORMAT[[1]] == 'AF'
	
	vcf$tech_rep_1[[1]] = vcf$tech_rep_1[[1]][vcf$FORMAT[[1]]]
	vcf$tech_rep_1 = as.character(vcf$tech_rep_1)
	vcf$tech_rep_1[1] = vcf$tech_rep_1[1] %>% stringr::str_split(',')
	vcf$tech_rep_1[1] = vcf$tech_rep_1[[1]][2] %>% as.numeric() 
	
	vcf$tech_rep_2[[1]] = vcf$tech_rep_2[[1]][vcf$FORMAT[[1]]]
	vcf$tech_rep_2 = as.character(vcf$tech_rep_2)
	vcf$tech_rep_2[1] = vcf$tech_rep_2[1] %>% stringr::str_split(',')
	vcf$tech_rep_2[1] = vcf$tech_rep_2[[1]][2] %>% as.numeric() 

	vcf$tech_rep_1 = as.numeric(vcf$tech_rep_1)
	vcf$tech_rep_2 = as.numeric(vcf$tech_rep_2)

	return(vcf)
}

vcf = purrr::map_dfr(.x=1:dim(vcf)[1], .f=get_AFs)
vcf$mean_MAF = ( vcf$tech_rep_1 + vcf$tech_rep_2 ) / 2
vcf$MAF_diff = abs( vcf$tech_rep_1 - vcf$tech_rep_2)

vcf = vcf %>% dplyr::select(-FORMAT)

# in the case of multi-amplicon hits, select the record with the lowest difference in reported MAFs between the two tech reps
vcf = vcf %>% dplyr::group_by(sample_id, vep_format) %>% dplyr::top_n(n=-1, wt=MAF_diff) %>% dplyr::ungroup()

print(vcf, width=Inf)
readr::write_tsv(vcf, snakemake@output[[1]])
