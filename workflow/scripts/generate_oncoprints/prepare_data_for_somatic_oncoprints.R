library(magrittr)
library(DBI)
library(RPostgres)

# TODO: move samples QC filtering upstream

prepare_data_for_oncoprint_generation = function (nonTP53_archival_variants, nonTP53_relapse_variants, TP53_variants, TP53_variant_clonality_status, gene_set_analysed, non_hgsoc_samples, samples_with_no_good_sequencing, samples_with_very_low_purity, output, analysis_type) {

	source('/Users/bradle02/.Renviron')
	source('functions.R')

	# read in tumour samples variants
	non_tp53_variants_archival = readr::read_tsv(nonTP53_archival_variants) %>% dplyr::mutate(type='archival')
	non_tp53_variants_relapse = readr::read_tsv(nonTP53_relapse_variants) %>% dplyr::mutate(type='relapse')
	
	# bind archival and relapse non-TP53 variants
	non_tp53_variants = rbind(non_tp53_variants_archival, non_tp53_variants_relapse)
	
	# remove suspected aritfacts
	non_tp53_variants = non_tp53_variants %>% dplyr::filter(vep_format!='chr17_41231352_T/C')

	# Set of annotated TP53 variants 
	tp53_variants = readr::read_tsv(TP53_variants) %>% 
			dplyr::mutate(SYMBOL='TP53')

	# restrict the TP53 variants analysed wrt patients based on the gene set analysis type
	if (gene_set_analysed == 'panel_6_28') {
	} else if (gene_set_analysed == 'panel_28_only') {
		panel_28_only_sequencing_metadata = readr::read_tsv('config/somatic_metadata_panel_28_only.tsv')
		tp53_variants = tp53_variants %>% dplyr::filter(fk_britroc_number %in% panel_28_only_sequencing_metadata$fk_britroc_number)
	}
	
	# read in TP53 variant clonality predictions
	tp53_variant_clonality = readr::read_tsv(TP53_variant_clonality_status)

	patients_with_no_clonality_calls = tp53_variant_clonality %>% dplyr::group_by(fk_britroc_number) %>% 
		dplyr::summarise(no_clonality_calls=all(is.na(classification))) %>% dplyr::filter(no_clonality_calls==TRUE) %>% 
		dplyr::pull(fk_britroc_number)

	join_TP53_variants_with_clonality_status = function(TP53_variant_set, TP53_variant_clonality_status, tumour_type) {

	## most samples have multiple TP53 mutations, so where applicable, we want to only select the clonal mutation

	# join the two TP53 variant tables
	TP53_variants_clonal_mutations = dplyr::inner_join(
		TP53_variant_set %>% dplyr::filter(type==tumour_type), 
		TP53_variant_clonality_status %>% # include variants without a classification and only those variants that appear in 'tumour_type' 
			dplyr::filter(classification=='clonal' | is.na(classification)) %>%
			dplyr::filter(stringr::str_interp('num_samples_with_variant_${tumour_type}')>0), 
		by=c('fk_britroc_number','#Uploaded_variation')
	) %>% dplyr::select(fk_britroc_number,Consequence,SYMBOL,type)
	
	return(TP53_variants_clonal_mutations)

	}
	
	tp53_variants_clonal_mutations_archival = join_TP53_variants_with_clonality_status(tp53_variants,tp53_variant_clonality,'archival')
	tp53_variants_clonal_mutations_relapse = join_TP53_variants_with_clonality_status(tp53_variants,tp53_variant_clonality,'relapse')

	# bind clonal mutations of different tumour types together
	tp53_variants_clonal_mutations = rbind(tp53_variants_clonal_mutations_archival, tp53_variants_clonal_mutations_relapse) %>% unique()

	# retain non-clonal variants for patients in which a clonal TP53 hasn't been found
	tp53_variants_no_clonal_mutation_archival = tp53_variants %>%
		dplyr::filter(type=='archival') %>% 
		dplyr::filter(!fk_britroc_number %in% tp53_variants_clonal_mutations_archival$fk_britroc_number) %>%
		dplyr::select(fk_britroc_number,Consequence,SYMBOL,type)
	
	tp53_variants_no_clonal_mutation_relapse = tp53_variants %>% 
		dplyr::filter(type=='relapse') %>%
		dplyr::filter(!fk_britroc_number %in% tp53_variants_clonal_mutations_relapse$fk_britroc_number) %>%
		dplyr::select(fk_britroc_number,Consequence,SYMBOL,type)
	
	tp53_variants_no_clonal_mutation = rbind(tp53_variants_no_clonal_mutation_archival, tp53_variants_no_clonal_mutation_relapse)

	# bind TP53 variants together
	tp53_variants = rbind(tp53_variants_clonal_mutations, tp53_variants_no_clonal_mutation)
	
	# reformat variant information
	non_tp53_variants = non_tp53_variants %>% dplyr::select(patient_id,Consequence,SYMBOL,type)
	tp53_variants = tp53_variants %>% dplyr::rename(patient_id=fk_britroc_number)

	# filter variants by type and by gene symbol on the basis of TMBP assigned pathogenicity status
	non_tp53_variants$SYMBOL %>% unique() %>% print()
	non_tp53_variants = non_tp53_variants %>% dplyr::filter(!SYMBOL %in% c('IFIT1P1'))

	if (gene_set_analysed=='panel_6_28') {

	# variants changed due to manual curation
	non_tp53_variants = non_tp53_variants %>% tibble::add_row(SYMBOL='BRCA2', patient_id=1, Consequence='frameshift_variant', type='archival')
	non_tp53_variants = non_tp53_variants %>% tibble::add_row(SYMBOL='BRCA2', patient_id=39, Consequence='frameshift_variant', type='archival')
	non_tp53_variants = non_tp53_variants %>% tibble::add_row(SYMBOL='TP53', patient_id=102, Consequence='frameshift_variant', type='archival')
	non_tp53_variants = non_tp53_variants %>% tibble::add_row(SYMBOL='BRCA2', patient_id=103, Consequence='frameshift_variant', type='archival')
	non_tp53_variants = non_tp53_variants %>% tibble::add_row(SYMBOL='BRCA2', patient_id=141, Consequence='frameshift_variant', type='archival')

	# variants removed due to manual curation
	non_tp53_variants = non_tp53_variants %>% dplyr::filter(!(SYMBOL == 'BRCA2' & patient_id==77 & type=='relapse'))
	non_tp53_variants = non_tp53_variants %>% dplyr::filter(!(SYMBOL == 'TP53' & patient_id==102 & type=='archival' & Consequence=='missense_variant')) # fixation artifact
	non_tp53_variants = non_tp53_variants %>% dplyr::filter(!(SYMBOL == 'BRCA2' & patient_id==112 & type=='archival'))
	non_tp53_variants = non_tp53_variants %>% dplyr::filter(!(SYMBOL == 'BRCA2' & patient_id==114 & type=='relapse'))

	} else if (gene_set_analysed == 'panel_28_only') {

	}

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
		
	relevant_samples = remove_non_relevant_samples(non_hgsoc_samples, samples_with_no_good_sequencing, samples_with_very_low_purity, snakemake@input['DNA_sample_file_path'], snakemake@input['slx_library_file_path'], snakemake@input['experiments_file_path'], analysis_type)

	# only retain variants in patients with at least one germline, relapse and archival sample
	all_variants = all_variants %>% dplyr::filter(patient_id %in% relevant_samples$fk_britroc_number)

	# reformat
	all_variants =
	  all_variants %>% dplyr::select(patient_id,SYMBOL,Consequence,type) %>%
	  unique()

	# remove synonymous and intron variants
	all_variants = 
	  all_variants %>% dplyr::filter(!Consequence %in% c('intron_variant','synonymous_variant')) %>%
	  unique()
	
	# this information is needed so that the oncoprint shows information even for patients without a mutation in any one gene
	# TODO: consider removing - purpose is not clear
	somatic_samples_with_no_mutations = dplyr::anti_join(
	  relevant_samples, all_variants, by=c('fk_britroc_number'='patient_id')
	)

	# restrict the TP53 variants analysed wrt patients based on the gene set analysis type
	if (gene_set_analysed == 'panel_6_28') {
	} else if (gene_set_analysed == 'panel_28_only') {
		panel_28_only_sequencing_metadata = readr::read_tsv('config/somatic_metadata_panel_28_only.tsv')
		somatic_samples_with_no_mutations = somatic_samples_with_no_mutations %>% dplyr::filter(fk_britroc_number %in% panel_28_only_sequencing_metadata$fk_britroc_number)
	}

	# reformat
	somatic_samples_with_no_mutations = somatic_samples_with_no_mutations %>%
	  dplyr::rename(patient_id='fk_britroc_number') %>%
	  dplyr::select(patient_id) %>% 
	  unique() %>% 
	  dplyr::mutate(variant_type='dummy', gene_symbol='PTEN', type='foo') %>% 
	  dplyr::select(patient_id,variant_type,gene_symbol,type)

	somatic_samples_with_mutations = all_variants %>% dplyr::select(patient_id,SYMBOL,Consequence,type)
	somatic_samples_with_no_mutations =  somatic_samples_with_no_mutations %>% dplyr::ungroup()
	colnames(somatic_samples_with_mutations) = c('patient_id','gene_symbol','variant_type','tumour_type')

	# bind data together
	somatic_samples_with_no_mutations = somatic_samples_with_no_mutations %>% dplyr::rename(
		'tumour_type'='type'
	)
	somatic_variants = rbind(somatic_samples_with_mutations, somatic_samples_with_no_mutations)

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
	  dplyr::group_by(patient_id,gene_symbol,tumour_type) %>%
	  dplyr::filter(dplyr::row_number()==1) %>% dplyr::ungroup()

	# reformat
	somatic_variants = somatic_variants %>% 
		tidyr::unite(gene_symbol_tumour_type, c(gene_symbol,tumour_type), remove=TRUE, sep='_') %>%
		dplyr::arrange(patient_id)

	readr::write_tsv(somatic_variants, output)					
}

prepare_data_for_oncoprint_generation (
	nonTP53_archival_variants=snakemake@input[['filtered_matched_and_paired_variants_archival']],
	nonTP53_relapse_variants=snakemake@input[['filtered_matched_and_paired_variants_relapse']],
	TP53_variants=snakemake@input[['filtered_TP53_variants_with_MAFs']],
	TP53_variant_clonality_status=snakemake@input[['clonality_status_of_TP53_variants']],
	gene_set_analysed=snakemake@wildcards$analysis_type,
	non_hgsoc_samples = snakemake@config[['non_hgsoc_samples']],
	samples_with_no_good_sequencing = snakemake@config[['samples_with_no_good_sequencing']],
	samples_with_very_low_purity = snakemake@config[['samples_with_very_low_purity']],
	output = snakemake@output[['data_for_somatic_oncoprint']],
	analysis_type='somatic'
)

#prepare_data_for_oncoprint_generation (
#	nonTP53_archival_variants='results/filtered_archival_vep_calls_octopus_joined.tsv',
#	nonTP53_relapse_variants='results/filtered_relapse_vep_calls_octopus_joined.tsv',
#	TP53_variants='results/final_tp53/filtered_TP53_variants_with_MAFs.tsv',
#	TP53_variant_clonality_status='results/final_tp53/TP53_variants_with_clonality_classifications.tsv',
#	non_hgsoc_samples = c('JBLAB-4114','JBLAB-4916','IM_249','IM_250','IM_234','IM_235','IM_236','IM_237','JBLAB-4271','IM_420','IM_262','JBLAB-4922','JBLAB-4923','IM_303','IM_290','IM_43','IM_293','IM_307','IM_308','IM_309','IM_424','IM_302','IM_303','IM_304','IM_305','JBLAB-19320','IM_61','IM_62','IM_63','IM_397','IM_302','IM_98','JBLAB-4210','IM_147','JBLAB-4216','IM_44'),
#	gene_set_analysed='HRD',
#	samples_with_no_good_sequencing =  c('IM_144','IM_435','IM_436','IM_158','IM_296','IM_373','IM_154','IM_297','IM_365','IM_432','IM_429','IM_368','IM_441'),
#	samples_with_very_low_purity = c('IM_1','IM_2','IM_3','IM_4','IM_20','IM_26','IM_27','IM_69','IM_86','IM_90','IM_93','IM_94','IM_173','IM_177','IM_179','IM_200','IM_241','IM_242','IM_417','IM_418','IM_419','IM_420','IM_221','IM_264','IM_329','IM_289','IM_308','IM_309','IM_338','IM_339','IM_340','IM_341','IM_342','IM_432','IM_372','IM_272','IM_392'),
#	output='somatic_variants_for_oncoprint.tsv'	
#)
