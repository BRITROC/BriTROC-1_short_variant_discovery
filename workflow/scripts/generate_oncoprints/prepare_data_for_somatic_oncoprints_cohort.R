library(magrittr)
library(DBI)
library(RPostgres)

prepare_data_for_oncoprint_generation = function (nonTP53_variants, TP53_variants, TP53_variant_clonality_status, non_hgsoc_samples, samples_with_no_good_sequencing, samples_with_very_low_purity, output, analysis_type) {

	source('/Users/bradle02/.Renviron')
	source('functions.R')

	print('goo1')

	britroc_con = make_connection_to_postgres_server('britroc1', 'jblab-db.cri.camres.org', 5432)
	clarity_con = make_connection_to_postgres_server('clarity', 'jblab-db.cri.camres.org', 5432)

	# read in tumour samples variants
	non_tp53_variants = readr::read_tsv(nonTP53_variants) %>% dplyr::filter(SYMBOL!='TP53')
	non_tp53_variants$HGVS_united = dplyr::coalesce(non_tp53_variants$HGVSp,non_tp53_variants$HGVSc)
	print(non_tp53_variants %>% dplyr::select(SYMBOL,HGVSc,HGVSp,POS), n=Inf, width=Inf)
	non_tp53_variants$HGVS_united = non_tp53_variants$HGVS_united %>% stringr::str_remove('ENS[A-Z0-9]+\\.[0-9]+:')

	non_tp53_variants = non_tp53_variants %>% dplyr::filter(!grepl('dup',HGVSc))

	print('goo2')

	# remove suspected aritfacts
	# evidence of existence in germline, v. likely spurious, inflated frequency compared to much higher confidence matched variant analysis

	#non_tp53_variants = non_tp53_variants %>% dplyr::filter(`#Uploaded_variation`!='chr17_41231352_T/C')
	non_tp53_variants = non_tp53_variants %>% dplyr::filter(`#Uploaded_variation`!='chr17_41245587_T/-')
 	non_tp53_variants = non_tp53_variants %>% dplyr::filter(`#Uploaded_variation`!='chr13_32954023_A/-') 

	# Set of annotated TP53 variants 
	tp53_variants = readr::read_tsv(TP53_variants) %>% 
			dplyr::mutate(SYMBOL='TP53')

	## most samples have multiple TP53 mutations, so where applicable, we want to only select the clonal mutation

	# read in TP53 variant clonality predictions
	tp53_variant_clonality = readr::read_tsv(TP53_variant_clonality_status)

	patients_with_no_clonality_calls = tp53_variant_clonality %>% dplyr::group_by(fk_britroc_number) %>% 
		dplyr::summarise(no_clonality_calls=all(is.na(classification))) %>% dplyr::filter(no_clonality_calls==TRUE) %>% 
		dplyr::pull(fk_britroc_number)

	join_TP53_variants_with_clonality_status = function(TP53_variant_set, TP53_variant_clonality_status) {

	# join the two TP53 variant tables
	TP53_variants_clonal_mutations = dplyr::inner_join(
		TP53_variant_set, 
		TP53_variant_clonality_status %>% # include variants without a classification
			dplyr::filter(classification=='clonal' | is.na(classification)) %>%
			dplyr::filter(stringr::str_interp('num_samples_with_variant_total')>0), 
		by=c('fk_britroc_number','#Uploaded_variation')
	) %>% dplyr::select(fk_britroc_number,Consequence,SYMBOL)
	
	return(TP53_variants_clonal_mutations)
	}
	
	tp53_variants_clonal_mutations = join_TP53_variants_with_clonality_status(tp53_variants,tp53_variant_clonality)

	# retain non-clonal variants for patients in which a clonal TP53 hasn't been found
	tp53_variants_no_clonal_mutation = tp53_variants %>%
		dplyr::filter(!fk_britroc_number %in% tp53_variants_clonal_mutations$fk_britroc_number) %>%
		dplyr::select(fk_britroc_number,Consequence,SYMBOL)

	# bind TP53 variants together
	tp53_variants = rbind(tp53_variants_clonal_mutations, tp53_variants_no_clonal_mutation)

	print('foo1')
	
	# reformat variant information
	non_tp53_variants = non_tp53_variants %>% dplyr::select(patient_id,Consequence,SYMBOL)
	tp53_variants = tp53_variants %>% dplyr::rename(patient_id=fk_britroc_number)

	# bind TP53 and nonTP53 variants together
	all_variants = rbind(tp53_variants, non_tp53_variants) %>% unique()

	# factorise variant consequence column
	all_variants$Consequence = factor(all_variants$Consequence,
	                                       levels=c(
	                                         'frameshift_variant',
	                                         'inframe_deletion',
	                                         'inframe_insertion',
	                                         'stop_gained',
	                                         'stop_gained,frameshift_variant',
	                                         'stop_gained,splice_region_variant',
	                                         'missense_variant',
	                                         'splice_donor_variant',
	                                         'splice_acceptor_variant',
	                                         'splice_region_variant,intron_variant',
						 'splice_region_variant,synonymous_variant'
	                                         ))

	# rearrange column
	all_variants = all_variants %>% dplyr::arrange(patient_id,Consequence)

	## remove non-relevant samples
	# TODO: remove SLX-13716 from the database
		
	relevant_samples = remove_non_relevant_samples(non_hgsoc_samples, samples_with_no_good_sequencing, samples_with_very_low_purity, britroc_con, clarity_con, analysis_type)

	# only retain variants in patients with relevant sample types available
	all_variants = all_variants %>% dplyr::filter(patient_id %in% relevant_samples$fk_britroc_number)

	# reformat
	all_variants =
	  all_variants %>% dplyr::select(patient_id,SYMBOL,Consequence) %>%
	  unique()
	
	print('foo2')

	# this information is needed so that the oncoprint shows information even for patients without a mutation in any one gene
	# TODO: consider removing - purpose is not clear
	somatic_samples_with_no_mutations = dplyr::anti_join(
	  relevant_samples, all_variants, by=c('fk_britroc_number'='patient_id')
	)
	
	# reformat
	somatic_samples_with_no_mutations = somatic_samples_with_no_mutations %>%
	  dplyr::rename(patient_id='fk_britroc_number') %>%
	  dplyr::select(patient_id) %>% 
	  unique() %>% 
	  dplyr::mutate(variant_type='dummy', gene_symbol='SHOE') %>% 
	  dplyr::select(patient_id,variant_type,gene_symbol)

	somatic_samples_with_mutations = all_variants %>% dplyr::select(patient_id,SYMBOL,Consequence)
	somatic_samples_with_no_mutations =  somatic_samples_with_no_mutations %>% dplyr::ungroup()
	colnames(somatic_samples_with_mutations) = c('patient_id','gene_symbol','variant_type')

	# bind data together
	somatic_variants = rbind(somatic_samples_with_mutations, somatic_samples_with_no_mutations)

	print('foo3')


	# add in high confidence germline variants
	# this is because octopus erroneously removes germline variants in most cases
	germline_variants = readr::read_tsv(snakemake@input[['germline_variants']])

	print(somatic_variants$patient_id %>% unique())
	somatic_variants$patient_id %>% unique() %>% length() %>% print()
	
	germline_variants = germline_variants %>% dplyr::filter(final_call_set==TRUE)
	germline_variants = germline_variants %>% dplyr::filter(fk_britroc_number %in% somatic_variants$patient_id)
	germline_variants = germline_variants %>% dplyr::filter(!is.na(variant_type))
	germline_variants = germline_variants %>% dplyr::rename(patient_id=fk_britroc_number)


	germline_variants = germline_variants %>% dplyr::select(patient_id,gene_symbol,variant_type) %>% unique()

	#print(somatic_variants)
	#print(germline_variants)

	somatic_variants = rbind(somatic_variants,germline_variants) %>% unique()
	#somatic_variants = somatic_variants %>% dplyr::filter(!gene_symbol %in% c('RAD51L3-RFFL','CTD-2196E14.3','SHOE'))

	somatic_variants$patient_id %>% unique() %>% print()
        somatic_variants$patient_id %>% unique() %>% length() %>% print()

	# remap variant type categories
	somatic_variants$variant_type = dplyr::recode(
	  somatic_variants$variant_type,
	  'frameshift,splice_region'='frameshift',
	  'stop_gained,frameshift_variant'='frameshift',
	  'frameshift_variant'='frameshift',
	  'inframe_deletion'='inframe indel',
	  'inframe_insertion'='inframe indel',
	  'stop_gained'='stop gained',
	  'stop_gained,splice_region_variant'='stop gained',
	  'stop_gained,splice_region'='stop gained',
	  'splice_acceptor'='splice region SNV',
	  'splice_region,intron'='splice region SNV',
	  'splice_region,synonymous'='splice region SNV',
	  'splice_acceptor_variant'='splice region SNV',
	  'splice_donor_variant'='splice region SNV',
	  'nonsynonymous'='missense',
	  'missense_variant'='missense',
	  "splice_region_variant,intron_variant" = 'splice region SNV',
	  "splice_region_variant,synonymous_variant"='splice region SNV',
	  'splice_region_variant'='splice region SNV',
	  'dummy'='' # make the dummy row blank
	)
 	
	#for simplicity for now - only one mutation per patient-gene group
	somatic_variants = somatic_variants %>% tibble::as_tibble() %>% 
	  dplyr::group_by(patient_id,gene_symbol) %>%
	  dplyr::filter(dplyr::row_number()==1) %>% dplyr::ungroup()

	readr::write_tsv(somatic_variants, output)					
}

print('shoe1')

prepare_data_for_oncoprint_generation (
	nonTP53_variants=snakemake@input[['filtered_non_TP53_variants']],
	TP53_variants=snakemake@input[['filtered_TP53_variants_with_MAFs']],
	TP53_variant_clonality_status=snakemake@input[['clonality_status_of_TP53_variants']],
	non_hgsoc_samples = snakemake@config[['non_hgsoc_samples']],
	samples_with_no_good_sequencing = snakemake@config[['samples_with_no_good_sequencing']],
	samples_with_very_low_purity = snakemake@config[['samples_with_very_low_purity']],
	output = snakemake@output[['data_for_somatic_oncoprint']],
	analysis_type='cohort'
)
