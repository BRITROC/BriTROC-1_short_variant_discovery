make_connection_to_postgres_server = function(database_name, host_name, port_number) {
	# aim: make a connection to the jblab postgres server using the DBI R package
	# input:
	#	database_name: The name of the database to connect to

	db_connection = DBI::dbConnect(
			RPostgres::Postgres(),
			dbname= database_name,
			host= host_name,
			port = port_number,
			user = jblab_db_username,
			password = jblab_db_password
	)

	return (db_connection)
}

remove_non_relevant_samples = function (non_hgsoc_samples, samples_with_no_good_sequencing, samples_with_very_low_purity, britroc_con, clarity_con, analysis_type) {
	
		# read in DNA sample and library information
		samples = dbReadTable(britroc_con, 'sample')
		
		slx_library = dbReadTable(britroc_con, 'slx_library') %>%
		       dplyr::filter(fk_slx!='SLX-13716') %>% # no data for this SLX
		       dplyr::filter(grepl('AA',fk_experiment)) # select for tam-seq experiments only
		
		slx_clarity = dbReadTable(clarity_con, 'slx')
		
		slx_library = dplyr::semi_join(slx_library, slx_clarity, by=c('fk_slx'='name')) # ensure slx is actually in clarity
		
		experiments = dbReadTable(britroc_con, 'experiment') %>% dplyr::filter(fk_amplicon_panel %in% c(1,6,10,28)) # only allow panel 6 or panel 28 amplicon panels
		slx_library = slx_library %>% dplyr::semi_join(experiments, by=c('fk_experiment'='name'))
		
		# retrieve analysis type
		sequenced_samples = dplyr::semi_join(samples, slx_library, by=c('name'='fk_sample'))
	
		# bind the poor samples together
		bad_samples = c(non_hgsoc_samples, samples_with_no_good_sequencing, samples_with_very_low_purity)
		
		# remove samples from downstream analyses
		sequenced_samples = sequenced_samples %>% dplyr::filter(!name %in% bad_samples)

		if (analysis_type=='somatic') {	
		# group sequenced samples by patient, and only retain patients with at least one germline, relapse and archival sample
		relevant_samples = sequenced_samples %>%
		       dplyr::group_by(fk_britroc_number,type) %>%
		       dplyr::summarise(n=dplyr::n()) %>%
		       tidyr::pivot_wider(names_from=type, values_from=n) %>%
		       dplyr::filter(!is.na(archival) & !is.na(relapse) & !is.na(germline)) # https://stackoverflow.com/questions/27197617/filter-data-frame-by-character-column-name-in-dplyr
		} else if (analysis_type=='germline') {
		relevant_samples = sequenced_samples %>%
		       dplyr::group_by(fk_britroc_number,type) %>%
		       dplyr::summarise(n=dplyr::n()) %>%
		       tidyr::pivot_wider(names_from=type, values_from=n) %>%
		       dplyr::filter(!is.na(germline))
		} else if (analysis_type=='cohort') {
		relevant_samples = sequenced_samples
		}

	return(relevant_samples)
		}

# TODO: formalise the process of mapping libraries to samples
# TODO: non-deterministic behaviour of string splitting which sometimes passes and sometimes fails
identify_variants_with_tech_rep_mismatch_in_joined_vcf_table = function(square_vcf) {
	# square_vcf: A joined vcf table where each column represents one library
	# output - sample_genotype_table: For each variant in square_vcf, a genotype for each sample where multiple libraries map onto single samples

	# create space for creating a table of genotypes
	genotype_table = square_vcf

	# extract the predicted genotype for each variant
	for (lib in libraries) {
		genotype_table[[lib]] = genotype_table[[lib]] %>% stringr::str_extract(pattern='[01\\|]+')
	}

	# go from a library genotype table to a sample genotype table - first by removing all library fields
	sample_genotype_table = genotype_table %>% dplyr::select(-dplyr::any_of(libraries)) 

	#TODO: revise this. Decide whether to keep this or not. For now keep but neutralise
	check_genotypes = function(row_index, genotype_table) {
		# genotype_table = A table of variants and their predicted genotypes for a set of libraries

		# output: A genotype table which only includes variants where the genotypes match between duplicate libraries
	
		new_genotype_table = genotype_table[row_index,]
		if (new_genotype_table[1,1] == new_genotype_table[1,2]) {
			return(new_genotype_table[1,1])
		} else {
			return(NA)
		}
	}

	# convert filtered/non-PASSed samples to the '0|0' genotype
	# FT refers to the FT field which scores each library on whether they passed all filters for a given variant
	square_vcf_ft = square_vcf
	if (snakemake@params[['includes_germline_variants']]==FALSE) {
		for (lib in libraries) {
			square_vcf_ft[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[17]]

			#print(lib)
			#print(square_vcf_ft[[lib]])

			genotype_table[[lib]] = ifelse(square_vcf_ft[[lib]] != 'PASS', '0|0', genotype_table[[lib]])
		}
	} else if (snakemake@params[['includes_germline_variants']]==TRUE) {
		square_vcf_ft_somatic = square_vcf_ft
		square_vcf_ft_germline_with_sh = square_vcf_ft
		square_vcf_ft_germline_without_sh = square_vcf_ft

		for (lib in libraries) {
			# without a somatic haplotype at that locus
			square_vcf_ft_germline_without_sh[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[13]]
			if (grepl('SOMATIC', square_vcf_ft[['INFO']]) %>% any()) { # test if somatic variants are present
				square_vcf_ft_somatic[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[17]]
				square_vcf_ft_germline_with_sh[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[14]]

				square_vcf_ft[[lib]] = dplyr::coalesce(square_vcf_ft_somatic[[lib]],square_vcf_ft_germline_with_sh[[lib]])
				square_vcf_ft[[lib]] = dplyr::coalesce(square_vcf_ft[[lib]],square_vcf_ft_germline_without_sh[[lib]])
			 # test if any germline variants occur on somatic haplotypes
			} else if ( (grepl('HSS', square_vcf_ft[['FORMAT']]) %>% any()) ) {
				# with somatic haplotype(s) at that locus
				square_vcf_ft_germline_with_sh[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[14]]
				
				square_vcf_ft[[lib]] = dplyr::coalesce(square_vcf_ft_germline_with_sh[[lib]],square_vcf_ft_germline_without_sh[[lib]])
			} else {
				square_vcf_ft[[lib]] = square_vcf_ft_germline_without_sh[[lib]]
			}
		genotype_table[[lib]] = ifelse(square_vcf_ft[[lib]] != 'PASS', '0|0', genotype_table[[lib]])
		}
	}

	# test that the two technical replicates are the same according to their predicted genotype
	for (sample in samples) {
		# subset to libraries relevant to a particular sample identifier
		relevant_libraries = genotype_table %>% dplyr::select(starts_with(sample)) 

		relevant_libraries$sample_genotype = purrr::map(.x=1:dim(relevant_libraries)[1], .f=check_genotypes, relevant_libraries) %>% unlist()

		# map the result of check_genotype onto a new sample-based table
		sample_genotype_table[[sample]] = relevant_libraries$sample_genotype
	}

	# filter for row records with at least one variant sample

	sample_genotype_table %>% dplyr::filter(POS==68301921) %>% print()
	sample_genotype_table = sample_genotype_table %>% dplyr::filter(dplyr::if_any(dplyr::all_of(samples), `%notin%`, c(NA,"0|0")))
	sample_genotype_table %>% dplyr::filter(POS==68301921) %>% print()

	sample_genotype_table = sample_genotype_table %>% dplyr::inner_join(square_vcf %>% dplyr::select(CHROM,POS,REF,ALT), by=c('CHROM','POS','REF','ALT'))

	return(sample_genotype_table)
}

implement_substitution_type_specific_filters = function(square_vcf) {
	# iterate over libraries and extract MAF and read depths
	if (snakemake@params[['includes_germline_variants']]==FALSE) {
	print('foo1')
		for (lib in libraries) {
			square_vcf_MAF[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[12]] %>%
				stringr::str_split(pattern=',') %>% data.table::transpose() %>% .[[2]] %>% stringr::str_replace(pattern='^.$', replacement='0.000') %>%
				as.numeric()
			square_vcf_depth[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[3]] %>%
				as.integer()
		}
	} else if (snakemake@params[['includes_germline_variants']]==TRUE) {
	print('foo2')
		square_vcf_MAF_somatic = square_vcf_MAF
		square_vcf_MAF_germline_with_sh = square_vcf_MAF # with somatic haplotype
		square_vcf_MAF_germline_without_sh = square_vcf_MAF # without somatic haplotype

		for (lib in libraries) {
			square_vcf_MAF_germline_with_sh[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[9]] %>%
				stringr::str_split(pattern=',') %>% data.table::transpose() %>% .[[2]] %>% stringr::str_replace(pattern='^.$', replacement='0.000') %>%
				as.numeric()
	
			square_vcf_MAF_germline_without_sh[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[8]] %>%
				stringr::str_split(pattern=',') %>% data.table::transpose() %>% .[[2]] %>% stringr::str_replace(pattern='^.$', replacement='0.000') %>%
				as.numeric()

			if (grepl('SOMATIC', square_vcf_MAF[['INFO']]) %>% any()) { # test if any of the variants are somatic
				print('foo3')
				square_vcf_MAF_somatic[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[12]] %>%
				stringr::str_split(pattern=',') %>% data.table::transpose() %>% .[[2]] %>% stringr::str_replace(pattern='^.$', replacement='0.000') %>%
					as.numeric()
				
				square_vcf_MAF[[lib]] = dplyr::coalesce(square_vcf_MAF_somatic[[lib]], square_vcf_MAF_germline_with_sh[[lib]])
				square_vcf_MAF[[lib]] = dplyr::coalesce(square_vcf_MAF[[lib]], square_vcf_MAF_germline_without_sh[[lib]])

			} else {
				print('foo4')
				square_vcf_MAF[[lib]] = dplyr::coalesce(square_vcf_MAF_germline_with_sh[[lib]], square_vcf_MAF_germline_without_sh[[lib]])
			}
			square_vcf_depth[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[3]] %>%
				as.integer()
		}
	}

	sample_weighted_MAF = square_vcf %>% dplyr::select(-all_of(libraries))
	
	for (sample in samples) {
		# extract the relevant depth and the relevant MAF for this sample
		relevant_depths = square_vcf_depth %>% dplyr::select(starts_with(sample))
		relevant_MAFs = square_vcf_MAF %>% dplyr::select(starts_with(sample))
	
		# weight MAFs by the respective depth at that locus 
		weighted_MAF_numerator = (relevant_depths[,1] * relevant_MAFs[,1]) + (relevant_depths[,2] * relevant_MAFs[,2])
		weighted_MAF_denominator = relevant_depths[,1] + relevant_depths[,2]
		weighted_MAF_numerator = weighted_MAF_numerator %>% dplyr::pull(sample)
		weighted_MAF_denominator = weighted_MAF_denominator %>% dplyr::pull(sample)
	
		# calculate the difference between MAFs of two technical replicates
		maf_diff = abs(relevant_MAFs[,1] - relevant_MAFs[,2]) %>% dplyr::pull(sample)
	
		# get the weighted MAF for each sample-variant combination
		sample_weighted_MAF[[sample]] = weighted_MAF_numerator / weighted_MAF_denominator
	
		# unite the REF and ALT columns
		sample_weighted_MAF = sample_weighted_MAF %>% tidyr::unite(col=substitution, REF,ALT, sep='>', remove=FALSE)
	
		# implement combined substitution type and MAF filters
		new_sample_genotype_column = ifelse(
			sample_weighted_MAF[[sample]] < snakemake@params[['C_to_G_maf_threshold']] & sample_weighted_MAF[['substitution']] %in% c('C>T','G>A'),
			'0|0', 
			sample_genotype_table[[sample]]
			)
		
		# implement substitution type and maf diff filters
		new_sample_genotype_column = ifelse(
			maf_diff >= snakemake@params[['C_to_G_maf_diff_threshold']] & sample_weighted_MAF[['substitution']] %in% c('G>A','C>T'),
			'0|0', 
			new_sample_genotype_column
			)
		
		sample_genotype_table[[sample]] = new_sample_genotype_column
	}

	# only retain variants which appear in at least one sample
	sample_genotype_table = sample_genotype_table %>% dplyr::filter(dplyr::if_any(all_of(samples), `%notin%`, c(NA,"0|0")))
	
	# use sample genotype table as a reference to filter other tables
	square_vcf_MAF = 
		dplyr::semi_join(
			square_vcf_MAF,
			sample_genotype_table,
			by=c('CHROM','POS','REF','ALT')
		)
	square_vcf_depth = 
		dplyr::semi_join(
			square_vcf_depth,
			sample_genotype_table,
			by=c('CHROM','POS','REF','ALT')
		)

	output_list = list(genotypes=sample_genotype_table,MAFs=square_vcf_MAF,depths=square_vcf_depth)

	return(output_list)
}
