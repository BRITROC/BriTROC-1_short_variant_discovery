# get tp53 table

library(magrittr)

archival_tp53_variants = readr::read_tsv('results/filtered_archival_vep_calls_octopus.tsv', na='NA')
relapse_tp53_variants = readr::read_tsv('results/filtered_relapse_vep_calls_octopus.tsv', na='NA')

all_tp53_variants = rbind(archival_tp53_variants,relapse_tp53_variants)

tp53_freq = all_tp53_variants %>% dplyr::group_by(fk_britroc_number,`#Uploaded_variation`) %>% dplyr::summarise(n=dplyr::n()) %>% dplyr::ungroup()

all_tp53_variants = dplyr::inner_join(all_tp53_variants,tp53_freq, by=c('fk_britroc_number','#Uploaded_variation'))

all_tp53_variants = all_tp53_variants %>% dplyr::arrange(fk_britroc_number,-n,type,sample_id)


MAFs = readr::read_tsv('results/tp53_collated_MAFs.tsv')

all_tp53_variants = dplyr::left_join(all_tp53_variants,MAFs, by=c('sample_id','#Uploaded_variation'='vep_format'))

# introduct changes relating to sample mislabelling

all_tp53_variants = all_tp53_variants %>% dplyr::filter(!sample_id %in% 
			c('JBLAB-4119','IM_92','IM_95','IM_116','IM_413')
		)

all_tp53_variants$fk_britroc_numer = ifelse(all_tp53_variants$sample_id=='IM_42','230',all_tp53_variants$fk_britroc_number)
all_tp53_variants$fk_britroc_numer = ifelse(all_tp53_variants$sample_id=='JBLAB-4968','233',all_tp53_variants$fk_britroc_number)

all_tp53_variants = all_tp53_variants %>% dplyr::group_by(sample_id,`#Uploaded_variation`) %>% dplyr::top_n(n=1, wt=QUAL) %>% dplyr::ungroup()

#colnames(all_tp53_variants) %>% print()

all_tp53_variants %>% dplyr::filter(sample_id=='IM_146') %>% dplyr::filter(`#Uploaded_variation`=='chr17_7577548_C/T') %>% 
	dplyr::select(sample_id,`#Uploaded_variation`,mean_MAF,MAF_diff) %>% print()

readr::write_tsv(
	all_tp53_variants %>% 
		dplyr::select(
			fk_britroc_number,sample_id,`#Uploaded_variation`,Consequence,CLIN_SIG,DOMAINS,type,QUAL,mean_MAF,MAF_diff
			) %>% unique(),
	'foo.tsv'
	)

all_tp53_variants = all_tp53_variants %>% unique()

mutation_freq = all_tp53_variants %>% dplyr::group_by(fk_britroc_number,type,`#Uploaded_variation`) %>% 
	dplyr::summarise(num_samples_with_variant=dplyr::n())

sample_freq = all_tp53_variants %>% dplyr::select(fk_britroc_number,type,sample_id) %>% unique() %>% dplyr::group_by(fk_britroc_number,type) %>%
	dplyr::summarise(num_samples=dplyr::n())

patient_level_table = dplyr::full_join(mutation_freq,sample_freq, by=c('fk_britroc_number','type'))

patient_level_table = patient_level_table %>%
	tidyr::pivot_wider(names_from=type, values_from=c('num_samples_with_variant','num_samples') )

patient_level_table = patient_level_table %>% dplyr::select(
	fk_britroc_number, `#Uploaded_variation`,num_samples_with_variant_archival,num_samples_archival, num_samples_with_variant_relapse,num_samples_relapse
	)

patient_level_table$num_samples_with_variant_archival = patient_level_table$num_samples_with_variant_archival %>% tidyr::replace_na(0)
patient_level_table$num_samples_with_variant_relapse = patient_level_table$num_samples_with_variant_relapse %>% tidyr::replace_na(0)

# fill in NAs
fill_in_nas = function(patient_id) {

	x = patient_level_table %>% dplyr::filter(fk_britroc_number==patient_id)

	if (x$num_samples_archival %>% is.na() %>% all == TRUE) {
		x$num_samples_archival = 0
	} else if (x$num_samples_archival %>% is.na() %>% any == TRUE) {
		num_samples_archival = x$num_samples_archival[!is.na(x$num_samples_archival)] %>% unique()
		x$num_samples_archival = num_samples_archival	
	}

	if (x$num_samples_relapse %>% is.na() %>% all == TRUE) {
		x$num_samples_relapse = 0
	} else if (x$num_samples_relapse %>% is.na() %>% any == TRUE) {
		num_samples_relapse = x$num_samples_relapse[!is.na(x$num_samples_relapse)] %>% unique()
		x$num_samples_relapse = num_samples_relapse	
	}

	return(x)
}

patient_level_table = purrr::map_dfr(patient_level_table$fk_britroc_number %>% unique(), .f=fill_in_nas)

patient_level_table = patient_level_table %>% dplyr::filter(!is.na(`#Uploaded_variation`))

patient_level_table$num_samples_with_variant_total = patient_level_table$num_samples_with_variant_archival + patient_level_table$num_samples_with_variant_relapse

patient_level_table = patient_level_table %>% dplyr::arrange(fk_britroc_number,-num_samples_with_variant_total)

#### Add MAF information

average_MAF_for_patient = all_tp53_variants %>% 
	dplyr::group_by(fk_britroc_number,`#Uploaded_variation`) %>%
	dplyr::summarise(mean_patient_MAF = mean(mean_MAF), mean_QUAL = mean(QUAL)) %>%
	dplyr::ungroup()

print(average_MAF_for_patient)

patient_level_table = patient_level_table %>% 
	dplyr::left_join(average_MAF_for_patient, by=c('fk_britroc_number','#Uploaded_variation')) %>%
	dplyr::arrange(fk_britroc_number,-mean_patient_MAF)
	
print(patient_level_table, width=Inf)

####

classify_clonality_status = function(patient_id) {

	x = patient_level_table %>% dplyr::filter(fk_britroc_number==patient_id)
	x$classification=NA

        # create ranks
	x$rank_num_samples = NA
	x$rank_num_samples[order(x$num_samples_with_variant_total, decreasing=TRUE)] = 1:nrow(x)

	x$rank_MAF = NA
	x$rank_MAF[order(x$mean_patient_MAF, decreasing=TRUE)] = 1:nrow(x)

	x$rank_QUAL = NA
	x$rank_QUAL[order(x$mean_QUAL, decreasing=TRUE)] = 1:nrow(x)	

	x$rank_total = x$rank_num_samples + x$rank_MAF + x$rank_QUAL

	# order by lowest summed rank
	x = x %>% dplyr::arrange(rank_total)
	
	#print(x %>% dplyr::select(num_samples_with_variant_total,mean_patient_MAF,mean_QUAL,rank_num_samples,rank_MAF,rank_QUAL, rank_total), width=Inf)
	#quit()

	if(x$mean_patient_MAF[1] < 0.07) {
		# if clonal mutation is less than 7% than do not attempt to make a call
	}
	else if(length(x$num_samples_with_variant_total) == 1) {
		x$classification[1] = 'clonal'
	}
	else {
	x$classification[1] = 'clonal'
	x$classification[2:length(x$classification)] = 'subclonal'
	}
	return(x)
}

patient_level_table = purrr::map_dfr(patient_level_table$fk_britroc_number %>% unique(), .f=classify_clonality_status)

patient_level_table = patient_level_table %>% dplyr::select(
	-rank_num_samples,-rank_MAF,-rank_QUAL,-rank_total		
)

#patient_level_table$num_samples_with_variant_total = NULL

readr::write_tsv(patient_level_table, 'foo2.tsv')
