
library(magrittr)

joined_vcf = readr::read_tsv('results/tumour_sample_vcfs_octopus/4.filtered3.vcf', comment='##')

joined_vcf = joined_vcf %>% tidyr::pivot_longer(
	cols=c(dplyr::starts_with('JBLAB'), dplyr::starts_with('IM')),
	names_to='library_id',
	values_to='format_values'
)

print(joined_vcf)

joined_vcf = joined_vcf %>% dplyr::mutate(sample_id = library_id %>% stringr::str_remove('_d$'))

joined_vcf = joined_vcf %>% dplyr::group_by(`#CHROM`,POS,REF,ALT)

joined_vcf = joined_vcf %>% dplyr::mutate(variant_id = dplyr::cur_group_id())

# split into two separate table
joined_vcf_samples = joined_vcf %>% dplyr::select(`#CHROM`,POS,REF,ALT,sample_id,INFO,variant_id) %>% unique() %>% dplyr::arrange(variant_id)
joined_vcf_libraries = joined_vcf %>% dplyr::select(`#CHROM`,POS,REF,ALT,library_id,sample_id,FORMAT,format_values,variant_id) %>% unique() %>% dplyr::arrange(variant_id)

# extract genotype information
joined_vcf_libraries = joined_vcf_libraries %>% dplyr::mutate(GT=format_values %>% stringr::str_extract(pattern='[0|1]+'))

# determine whether genotypes for libraries within a sample match

do_genotypes_match = function(sample,variant) {

	joined_vcf_libraries = joined_vcf_libraries %>% 
		dplyr::filter(variant_id==variant) %>% 
		dplyr::filter(sample_id==sample)

	num_genotypes = joined_vcf_libraries$GT %>% unique() %>% length() 

	if (num_genotypes==1) {
		return (joined_vcf_libraries$GT[1])
	} else {
		return(NA)
	}
}

print(joined_vcf_samples)
print(joined_vcf_libraries, width=Inf)

sample_genotypes = purrr::map2(.x=joined_vcf_samples$sample_id, .y=joined_vcf_samples$variant_id, .f=do_genotypes_match) %>% unlist()

joined_vcf_samples$genotype = sample_genotypes
joined_vcf_samples$non_wt_genotype = ifelse(joined_vcf_samples$genotype != '0|0', TRUE, FALSE)
joined_vcf_samples$non_wt_genotype = ifelse(is.na(joined_vcf_samples$genotype), FALSE, joined_vcf_samples$non_wt_genotype)

print(joined_vcf_samples %>% dplyr::select(-INFO), n=Inf)

# set up separate variant level table

variants = tibble::tibble(variant_id=joined_vcf_samples$variant_id %>% unique())

is_variant_valid_wrt_genotype = function(variant) {
	joined_vcf_samples = joined_vcf_samples %>% dplyr::filter(variant_id==variant)
	
	if (joined_vcf_samples$non_wt_genotype %>% any()) {
		return(TRUE)
	} else {
		return(FALSE)
	}	
}

variants$keep_variant = purrr::map(.x=joined_vcf_samples$variant_id %>% unique, .f=is_variant_valid_wrt_genotype) %>% unlist()
variants = variants %>% dplyr::filter(keep_variant==TRUE) %>% dplyr::select(-keep_variant)

print(variants, n=Inf)

joined_vcf_samples = dplyr::inner_join(joined_vcf_samples,variants,by='variant_id')
joined_vcf_libraries = dplyr::inner_join(joined_vcf_libraries,variants,by='variant_id')

print(joined_vcf_samples %>% dplyr::select(-INFO), width=Inf, n=Inf)
#print(joined_vcf_libraries)

joined_vcf_samples %>% dplyr::group_by(variant_id,non_wt_genotype) %>% dplyr::summarise(n=dplyr::n()) %>% dplyr::arrange(-non_wt_genotype,-n) %>% print(n=Inf)

#print(joined_vcf %>% dplyr::sample_n(size=1), width=Inf)
