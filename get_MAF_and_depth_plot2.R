
library(magrittr)
library(ggplot2)
library(gridExtra)
library(patchwork)

archival_variants = readr::read_tsv('tmp_annotations_joined_archival.tsv') %>% dplyr::select(patient_id,CHROM,POS,REF,ALT) %>% dplyr::arrange(patient_id)
metadata = readr::read_tsv('config/matched_somatic_metadata.tsv')

interval_list = readr::read_tsv('resources/intersected_panel_6_28_amplicons.targets.interval_list', comment='@', col_names=FALSE)
interval_list$gene = interval_list$X5 %>% stringr::str_extract('[A-Z0-9]+')
interval_list$length = abs(interval_list$X2 - interval_list$X3)

interval_list = interval_list %>% dplyr::group_by(gene) %>% dplyr::summarise(sum_length=sum(length))
 
print(interval_list)
#print(archival_variants)

`%notin%` = function(x,y) !(x %in% y)

pdf('maf_temp.pdf')

generate_plot = function(patient_id, CHROM, POS, REF, ALT) {

	#print(patient_id)
	#print(CHROM)
	#print(POS)
	#print(REF)
	#print(ALT)

	if (nchar(REF) > 1) {
		POS = POS - 1
	}

	library_MAFs = readr::read_tsv(
		stringr::str_interp('results/tumour_sample_vcfs_octopus/${patient_id}.library_MAFs.vcf')
		)

	library_depths = readr::read_tsv(
		stringr::str_interp('results/tumour_sample_vcfs_octopus/${patient_id}.library_depths.vcf')
		)

	# identify the relevant record
	library_MAFs = library_MAFs %>% dplyr::filter(CHROM==!!CHROM,POS==!!POS,REF==!!REF,ALT==!!ALT)
	library_depths = library_depths %>% dplyr::filter(CHROM==!!CHROM,POS==!!POS,REF==!!REF,ALT==!!ALT)

	# identify relevant samples
	archival_samples = metadata %>% dplyr::filter(fk_britroc_number==!!patient_id, type=='archival') %>% dplyr::pull(fk_sample) %>% unique()
	relapse_samples = metadata %>% dplyr::filter(fk_britroc_number==!!patient_id, type=='relapse') %>% dplyr::pull(fk_sample) %>% unique()

	# identify the sample(s) which contains the variant
	sample_genotypes = readr::read_tsv(
		stringr::str_interp('results/tumour_sample_vcfs_octopus/${patient_id}.sample_genotypes.vcf'),
	col_types = 
		readr::cols(
			POS='c',
			ALT='c'	
		)
	)
	#print('sample_genotypes')
	#print(sample_genotypes)
		# identify the relevant record
	sample_genotypes = sample_genotypes %>% dplyr::filter(CHROM==!!CHROM,POS==!!POS,REF==!!REF,ALT==!!ALT)

	#print(sample_genotypes, width=Inf)

	archival_samples_with_variant = sample_genotypes %>% 
			dplyr::select(all_of(archival_samples)) %>% 
			tidyr::pivot_longer(all_of(archival_samples), names_to='sample', values_to='genotype') %>%
			dplyr::filter(!genotype %in% c(NA,"0|0")) %>% 
			.$sample %>% unique()

	#print(archival_samples_with_variant)

	MAF_table = tibble::tibble()

	if (length(archival_samples_with_variant) != 0) {
		archival_libraries_with_variant = c(archival_samples_with_variant  ,paste(archival_samples_with_variant,'d',sep='_'))

		archival_libraries_with_variant_MAF = library_MAFs %>% dplyr::select(all_of(archival_libraries_with_variant)) %>%
		tidyr::pivot_longer(cols = all_of(archival_libraries_with_variant), names_to='library', values_to='MAF') %>%
		dplyr::mutate(lib_type='archival_libs_variant')

		archival_libraries_with_variant_depth = library_depths %>% dplyr::select(all_of(archival_libraries_with_variant)) %>%
		tidyr::pivot_longer(cols = all_of(archival_libraries_with_variant), names_to='library', values_to='depth')

		archival_libraries_with_variant_MAF = archival_libraries_with_variant_MAF %>% 
			dplyr::inner_join(archival_libraries_with_variant_depth, by='library')

		MAF_table = rbind(MAF_table,archival_libraries_with_variant_MAF)
		MAF_table$sample = MAF_table$library %>% stringr::str_remove('_d')
		MAF_table$library  = MAF_table$library %>% stringr::str_replace(pattern=archival_samples_with_variant, 'rep')
		MAF_table$patient_id = patient_id
		MAF_table$CHROM = CHROM
		MAF_table$POS = POS
		MAF_table$REF = REF
		MAF_table$ALT = ALT
		MAF_table = MAF_table %>% dplyr::select(-depth)

		print(MAF_table)

		return(MAF_table)
	} 
}

x = purrr::pmap_dfr(archival_variants %>% dplyr::select(patient_id, CHROM, POS, REF, ALT), generate_plot )

x = x %>% tidyr::pivot_wider(id_cols=c(patient_id, CHROM, POS, REF, ALT, sample), names_from=library, values_from=MAF)

x$MAF_diff = abs(x$rep - x$rep_d) 
x$mean_MAF = (x$rep + x$rep_d) /2 

x = x %>% dplyr::filter(MAF_diff < 0.15)

print(x)

archival_variants = readr::read_tsv('tmp_annotations_joined_archival.tsv')

x2 = x %>% dplyr::inner_join(archival_variants, by=c('patient_id','CHROM','POS','REF','ALT')) %>% dplyr::select('patient_id','CHROM','POS','REF','ALT', 'sample', 'SYMBOL','mean_MAF')

x2 = x2 %>% tidyr::unite(substitution, REF, ALT, remove=TRUE, sep='>')

x2$substitution = dplyr::recode(x2$substitution,
                                               'G>T'='C>A',
                                               'G>C'='C>G',
                                               'G>A'='C>T',
                                               'A>T'='T>A',
                                               'A>G'='T>C',
                                               'A>C'='T>G'
                                               )

print(x2)

print(x2$substitution %>% table())

ggplot(x2, aes(x=substitution)) + geom_bar(fill="#FF9999", colour="black") + 
	ggtitle('BritROC-1: Substitutions detected in archival samples') +
	labs(subitle='Mutation found in both technical replicates')

#ggsave('foo.png', device='png')

#quit()

print(x2$mean_MAF %>% summary())

ggplot(x2, aes(x=mean_MAF)) + geom_density() + ggtitle('BritROC-1: Substitutions detected in archival samples')

ggsave('foo.png', device='png')

quit()

x3 = x2 %>% dplyr::group_by(SYMBOL) %>% dplyr::summarise(n=dplyr::n())

print(interval_list)

x3 = x3 %>% dplyr::inner_join(interval_list, by=c('SYMBOL'='gene'))

print(x3)

x3$normalised_n = x3$n / as.integer(x3$sum_length)

ggplot(x3, aes(x=SYMBOL, y=normalised_n)) + geom_col()

ggsave('foo.png', device='png') 






#dplyr::anti_join(x, x2, by=c('patient_id','CHROM','POS','REF','ALT', 'sample'))

#ggsave(plot = marrangeGrob(p, nrow=1, ncol=1), 'maf_temp.pdf', device='pdf')
