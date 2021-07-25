
library(magrittr)
library(ggplot2)
library(gridExtra)
library(patchwork)

archival_variants = readr::read_tsv('tmp_annotations_joined_archival.tsv') %>% dplyr::select(patient_id,CHROM,POS,REF,ALT) %>% dplyr::arrange(patient_id)
metadata = readr::read_tsv('config/matched_somatic_metadata.tsv')

#print(archival_variants)

`%notin%` = function(x,y) !(x %in% y)

pdf('maf_temp_new.pdf')

generate_plot = function(patient_id, CHROM, POS, REF, ALT) {

	print(patient_id)
	print(CHROM)
	print(POS)
	print(REF)
	print(ALT)

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
	print('sample_genotypes')
	print(sample_genotypes)
		# identify the relevant record
	sample_genotypes = sample_genotypes %>% dplyr::filter(CHROM==!!CHROM,POS==!!POS,REF==!!REF,ALT==!!ALT)

	print(sample_genotypes)

	archival_samples_with_variant = sample_genotypes %>% 
			dplyr::select(all_of(archival_samples)) %>% 
			tidyr::pivot_longer(all_of(archival_samples), names_to='sample', values_to='genotype') %>%
			dplyr::filter(!genotype %in% c(NA,"0|0")) %>% 
			.$sample %>% unique()

	relapse_samples_with_variant = sample_genotypes %>% 
			dplyr::select(all_of(relapse_samples)) %>% 
			tidyr::pivot_longer(all_of(relapse_samples), names_to='sample', values_to='genotype') %>%
			dplyr::filter(!genotype %in% c(NA,"0|0")) %>% 
			.$sample %>% unique()
	

	
	#if(dim(archival_samples_with_variant)[1] != 0) {
	#	archival_samples_with_variant = colnames(archival_samples_with_variant)
	#} else {
	#	archival_samples_with_variant = NA
	#}

	archival_samples_without_variant = sample_genotypes %>% 
			dplyr::select(all_of(archival_samples)) %>% 
			tidyr::pivot_longer(all_of(archival_samples), names_to='sample', values_to='genotype') %>%
			dplyr::filter(genotype %in% c(NA,"0|0")) %>% 
			.$sample %>% unique()

	relapse_samples_without_variant = sample_genotypes %>% 
			dplyr::select(all_of(relapse_samples)) %>% 
			tidyr::pivot_longer(all_of(relapse_samples), names_to='sample', values_to='genotype') %>%
			dplyr::filter(genotype %in% c(NA,"0|0")) %>% 
			.$sample %>% unique()

	#archival_samples_without_variant = sample_genotypes %>% 
	#		dplyr::select(all_of(archival_samples)) %>% 
	#		dplyr::filter(dplyr::if_any(all_of(archival_samples), `%in%`, c(NA,"0|0")))

	#relapse_samples_without_variant = sample_genotypes %>% 
	#		dplyr::select(all_of(relapse_samples)) %>% 
	#		dplyr::filter(dplyr::if_any(all_of(relapse_samples), `%in%`, c(NA,"0|0")))

	print(archival_samples_with_variant)
	print(archival_samples_without_variant)
	print(relapse_samples_with_variant)
	print(relapse_samples_without_variant)

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
	} 

	if (length(relapse_samples_with_variant) != 0) {
		relapse_libraries_with_variant = c(relapse_samples_with_variant  ,paste(relapse_samples_with_variant,'d',sep='_'))

		relapse_libraries_with_variant_MAF = library_MAFs %>% dplyr::select(all_of(relapse_libraries_with_variant)) %>%
		tidyr::pivot_longer(cols = all_of(relapse_libraries_with_variant), names_to='library', values_to='MAF') %>%
		dplyr::mutate(lib_type='relapse_libs_variant')

		relapse_libraries_with_variant_depth = library_depths %>% dplyr::select(all_of(relapse_libraries_with_variant)) %>%
		tidyr::pivot_longer(cols = all_of(relapse_libraries_with_variant), names_to='library', values_to='depth')

		relapse_libraries_with_variant_MAF = relapse_libraries_with_variant_MAF %>% 
			dplyr::inner_join(relapse_libraries_with_variant_depth, by='library')

		MAF_table = rbind(MAF_table,relapse_libraries_with_variant_MAF)
	}

	if (length(archival_samples_without_variant) != 0) {
		archival_libraries_without_variant = c(archival_samples_without_variant  ,paste(archival_samples_without_variant,'d',sep='_'))

		archival_libraries_without_variant_MAF = library_MAFs %>% dplyr::select(all_of(archival_libraries_without_variant)) %>%
		tidyr::pivot_longer(cols = all_of(archival_libraries_without_variant), names_to='library', values_to='MAF') %>%
		dplyr::mutate(lib_type='archival_libs_no_variant')

		archival_libraries_without_variant_depth = library_depths %>% dplyr::select(all_of(archival_libraries_without_variant)) %>%
		tidyr::pivot_longer(cols = all_of(archival_libraries_without_variant), names_to='library', values_to='depth')

		archival_libraries_without_variant_MAF = archival_libraries_without_variant_MAF %>% 
			dplyr::inner_join(archival_libraries_without_variant_depth, by='library')

		MAF_table = rbind(MAF_table,archival_libraries_without_variant_MAF)
	}

	if (length(relapse_samples_without_variant) != 0) {
		relapse_libraries_without_variant = c(relapse_samples_without_variant  ,paste(relapse_samples_without_variant,'d',sep='_'))

		relapse_libraries_without_variant_MAF = library_MAFs %>% dplyr::select(all_of(relapse_libraries_without_variant)) %>%
		tidyr::pivot_longer(cols = all_of(relapse_libraries_without_variant), names_to='library', values_to='MAF') %>%
		dplyr::mutate(lib_type='relapse_libs_no_variant')

		relapse_libraries_without_variant_depth = library_depths %>% dplyr::select(all_of(relapse_libraries_without_variant)) %>%
		tidyr::pivot_longer(cols = all_of(relapse_libraries_without_variant), names_to='library', values_to='depth')

		relapse_libraries_without_variant_MAF = relapse_libraries_without_variant_MAF %>% 
			dplyr::inner_join(relapse_libraries_without_variant_depth, by='library')

		MAF_table = rbind(MAF_table,relapse_libraries_without_variant_MAF)
	}

	print(MAF_table)

	if (!'lib_type' %in% colnames(MAF_table)) {
		print('no lib type')
		return(NULL)
	}

	p_MAF = ggplot(MAF_table, aes(x=lib_type, y=MAF)) + geom_point(
			aes(colour=lib_type),
			position = position_dodge2(preserve = "single", width=0.2)
			) + ylim(c(0,1)) + ggtitle(paste('patient',patient_id, CHROM, POS, REF, ALT, sep='_'))

	p_depth = ggplot(MAF_table, aes(x=lib_type, y=depth)) + geom_point(
		aes(colour=lib_type),
		position = position_dodge2(preserve = "single", width=0.2)
		) + ggtitle(paste('patient',patient_id, CHROM, POS, REF, ALT, sep='_'))

	p_combined <- p_MAF / p_depth
	p_combined
	#return(p_combined)

	#ggsave('maf_temp.png', device='png')
}

purrr::pmap(archival_variants %>% dplyr::select(patient_id, CHROM, POS, REF, ALT), generate_plot )

dev.off()

#ggsave(plot = marrangeGrob(p, nrow=1, ncol=1), 'maf_temp.pdf', device='pdf')
