
library(magrittr)

square_vcf = readr::read_tsv('105.filtered3.vcf', comment='##')

# remove final column in the table which corresponds to the normal sample
square_vcf = square_vcf[,1:length(colnames(square_vcf))-1]

square_vcf = square_vcf %>% dplyr::filter(QUAL>500)
#square_vcf = square_vcf[10,]
genotype_table = square_vcf

libraries = genotype_table %>% dplyr::select(dplyr::starts_with('IM'), dplyr::starts_with('JBLAB')) %>% colnames()

for (lib in libraries) {
	genotype_table[[lib]] = genotype_table[[lib]] %>% stringr::str_extract(pattern='[01\\|]+')
}

# go from a library genotype table to a sample genotype table
sample_genotype_table = genotype_table %>% dplyr::select(-dplyr::starts_with('IM'), -dplyr::starts_with('JBLAB'))

samples = libraries %>% stringr::str_remove(pattern='_d') %>% unique()

check_genotypes = function(row_index, genotype_table) {
	
	new_genotype_table = genotype_table[row_index,]
	if (new_genotype_table[1,1] == new_genotype_table[1,2]) {
		return(new_genotype_table[1,1])
	} else {
		return(NA)
	}
}

# test that the two technical replicates are the same
for (sample in samples) {
	relevant_libraries = genotype_table %>% dplyr::select(starts_with(sample))
	#relevant_libraries$sample_genotype = relevant_libraries[,1] == relevant_libraries[,2] 

	relevant_libraries$sample_genotype = purrr::map(.x=1:dim(relevant_libraries)[1], .f=check_genotypes, relevant_libraries) %>% unlist()

	#relevant_libraries$sample_genotype = ifelse(relevant_libraries[,1] == relevant_libraries[,2], relevant_libraries %>% dplyr::select(ends_with('_d')), NA)
	sample_genotype_table[[sample]] = relevant_libraries$sample_genotype
}

# filter row records with at least one variant sample
`%notin%` = function(x,y) !(x %in% y)
sample_genotype_table = sample_genotype_table %>% dplyr::filter(dplyr::if_any(c(starts_with('IM'), starts_with('JBLAB')), `%notin%`, c(NA,"0|0")))

square_vcf = square_vcf %>% dplyr::inner_join(sample_genotype_table %>% dplyr::select('#CHROM',POS,REF,ALT), by=c('#CHROM','POS','REF','ALT'))
print(square_vcf, width=Inf)

# get MAFs for each sample

#square_vcf$IM_403 %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[3]] #%>% 
	#stringr::str_split(pattern=',') %>% data.table::transpose() %>% .[[2]] %>% stringr::str_replace(pattern='^.$', replacement='0.000')

square_vcf_MAF = square_vcf
square_vcf_depth = square_vcf

for (lib in libraries) {
	square_vcf_MAF[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[12]] %>%
		stringr::str_split(pattern=',') %>% data.table::transpose() %>% .[[2]] %>% stringr::str_replace(pattern='^.$', replacement='0.000') %>%
		as.numeric()
	square_vcf_depth[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[3]] %>%
		as.integer()
}

print(square_vcf_MAF)
print(square_vcf_depth)

sample_weighted_MAF = square_vcf %>% dplyr::select(-dplyr::starts_with('IM'), -dplyr::starts_with('JBLAB'))

for (sample in samples) {
	print(sample)
	relevant_depths = square_vcf_depth %>% dplyr::select(starts_with(sample))
	relevant_MAFs = square_vcf_MAF %>% dplyr::select(starts_with(sample))

	#relevant_libraries$sample_genotype = purrr::map(.x=1:dim(relevant_libraries)[1], .f=check_genotypes, relevant_libraries) %>% unlist()

	weighted_MAF_numerator = (relevant_depths[,1] * relevant_MAFs[,1]) + (relevant_depths[,2] * relevant_MAFs[,2])
	weighted_MAF_denominator = relevant_depths[,1] + relevant_depths[,2]
	weighted_MAF_numerator = weighted_MAF_numerator %>% dplyr::pull(sample)
	weighted_MAF_denominator = weighted_MAF_denominator %>% dplyr::pull(sample)

	sample_weighted_MAF[[sample]] = weighted_MAF_numerator / weighted_MAF_denominator
}

print(square_vcf_MAF %>% dplyr::select(starts_with('IM'), starts_with('JBLAB')), n=Inf)
print(square_vcf_depth %>% dplyr::select(starts_with('IM'), starts_with('JBLAB')), n=Inf)
print(sample_weighted_MAF %>% dplyr::select(starts_with('IM'), starts_with('JBLAB')), width=Inf, n=Inf)
print(sample_genotype_table %>% dplyr::select(-FORMAT), width=Inf, n=Inf)
