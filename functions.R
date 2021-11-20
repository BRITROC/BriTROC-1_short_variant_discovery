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
		
		experiments = dbReadTable(britroc_con, 'experiment') %>% dplyr::filter(fk_amplicon_panel %in% c(6,28)) # only allow panel 6 or panel 28 amplicon panels
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

