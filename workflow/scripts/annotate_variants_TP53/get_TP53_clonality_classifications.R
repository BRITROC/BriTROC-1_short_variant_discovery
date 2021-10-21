# classify TP53 variants
library(magrittr)

# read in TP53 variants with MAF information
tp53_variants_with_MAFs = readr::read_tsv(snakemake@input[['filtered_TP53_variants_with_MAFs']])
   
tp53_variants_with_MAFs = tp53_variants_with_MAFs %>% unique()

# derive the frequency of each variant for each patient and for each sample type
mutation_freq = tp53_variants_with_MAFs %>% 
		dplyr::group_by(fk_britroc_number,type,`#Uploaded_variation`) %>% 
		dplyr::summarise(num_samples_with_variant=dplyr::n())

# determine the number of tumour samples of each type for each patient
sample_freq = tp53_variants_with_MAFs %>% 
		dplyr::select(fk_britroc_number,type,sample_id) %>% 
		unique() %>% 
		dplyr::group_by(fk_britroc_number,type) %>%
		dplyr::summarise(num_samples=dplyr::n())

# join the previously generated tables 
patient_level_table = dplyr::full_join(mutation_freq,sample_freq, by=c('fk_britroc_number','type'))

# reformat table to make it wider 
patient_level_table = patient_level_table %>%
	tidyr::pivot_wider(names_from=type, values_from=c('num_samples_with_variant','num_samples') )

# select columns of interest
patient_level_table = patient_level_table %>% dplyr::select(
	fk_britroc_number, `#Uploaded_variation`,num_samples_with_variant_archival,num_samples_archival, num_samples_with_variant_relapse,num_samples_relapse
	)

# replace NA values with 0 for the num samples with variant columns
patient_level_table$num_samples_with_variant_archival = patient_level_table$num_samples_with_variant_archival %>% tidyr::replace_na(0)
patient_level_table$num_samples_with_variant_relapse = patient_level_table$num_samples_with_variant_relapse %>% tidyr::replace_na(0)

# fill in NAs
fill_in_nas = function(patient_id) {
	# replace nas for the num sample columns with the correct number of samples for that patient and tumour type

	# filter for patient
	x = patient_level_table %>% dplyr::filter(fk_britroc_number==patient_id)

	# replace na values for num samples 
	if (x$num_samples_archival %>% is.na() %>% all == TRUE) { # if all variants are NA for this tumour type, set num_samples to 0
		x$num_samples_archival = 0
	} else if (x$num_samples_archival %>% is.na() %>% any == TRUE) { # if some are NA, replace NA value with num samples for this sample type
		num_samples_archival = x$num_samples_archival[!is.na(x$num_samples_archival)] %>% unique()
		x$num_samples_archival = num_samples_archival	
	}

	# repeat for all tumour types
	if (x$num_samples_relapse %>% is.na() %>% all == TRUE) {
		x$num_samples_relapse = 0
	} else if (x$num_samples_relapse %>% is.na() %>% any == TRUE) {
		num_samples_relapse = x$num_samples_relapse[!is.na(x$num_samples_relapse)] %>% unique()
		x$num_samples_relapse = num_samples_relapse	
	}

	return(x)
}

# replace NAs for the number of samples for given sample type and patient
patient_level_table = purrr::map_dfr(patient_level_table$fk_britroc_number %>% unique(), .f=fill_in_nas)

# remove empty records
patient_level_table = patient_level_table %>% dplyr::filter(!is.na(`#Uploaded_variation`))

# determine the total number of samples with variant for each patient
patient_level_table$num_samples_with_variant_total = patient_level_table$num_samples_with_variant_archival + patient_level_table$num_samples_with_variant_relapse

# arrange the table according to patient and number of samples with the variant for each variant
patient_level_table = patient_level_table %>% dplyr::arrange(fk_britroc_number,-num_samples_with_variant_total)

#### Add MAF information

# Calculate summary MAF statistics
average_MAF_for_patient = tp53_variants_with_MAFs %>% 
	dplyr::group_by(fk_britroc_number,`#Uploaded_variation`) %>%
	dplyr::summarise(mean_patient_MAF = mean(mean_MAF), mean_QUAL = mean(QUAL)) %>%
	dplyr::ungroup()

# join MAF information into the main table
patient_level_table = patient_level_table %>% 
	dplyr::left_join(average_MAF_for_patient, by=c('fk_britroc_number','#Uploaded_variation')) %>%
	dplyr::arrange(fk_britroc_number,-mean_patient_MAF)
	
# classify clonality status of each TP53 mutation on a patient-by-patient basis
classify_clonality_status = function(patient_id) {

	# filter for patient ID
	x = patient_level_table %>% 
		dplyr::filter(fk_britroc_number==patient_id)

	# generate new classification and rank sum column
	x$classification=NA
	x$rank_num_samples = NA

	# rank variants by total number of samples in which they appear
	x$rank_num_samples[order(x$num_samples_with_variant_total, decreasing=TRUE)] = 1:nrow(x)

	# rank variants by their mutant allele fraction
	x$rank_MAF = NA
	x$rank_MAF[order(x$mean_patient_MAF, decreasing=TRUE)] = 1:nrow(x)

	# rank variants by their quality score
	x$rank_QUAL = NA
	x$rank_QUAL[order(x$mean_QUAL, decreasing=TRUE)] = 1:nrow(x)	

	# add a penalty for variants which don't appear in all tumour sample types
	x$tumour_type_penalty = NA

	x$tumour_type_penalty = ifelse(x$num_samples_with_variant_archival > 0 & x$num_samples_with_variant_relapse > 0 , 0, 4) # a penalty of 4
	x$tumour_type_penalty = ifelse(x$num_samples_archival == 0 | x$num_samples_relapse == 0 , 0, x$tumour_type_penalty)

	# calculate the sum of ranks
	x$rank_total = x$rank_num_samples + x$rank_MAF + x$rank_QUAL + x$tumour_type_penalty

	# order by lowest summed rank
	x = x %>% dplyr::arrange(rank_total)

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

# make a clonality classification on each TP53 mutation 
patient_level_table = purrr::map_dfr(patient_level_table$fk_britroc_number %>% unique(), .f=classify_clonality_status)

# remove redundant columns from the table before saving
patient_level_table = patient_level_table %>% dplyr::select(
	-rank_num_samples,-rank_MAF,-rank_QUAL,-rank_total		
)

# write table to file
readr::write_tsv(patient_level_table, snakemake@output[['TP53_variants_classified_by_clonality']])
