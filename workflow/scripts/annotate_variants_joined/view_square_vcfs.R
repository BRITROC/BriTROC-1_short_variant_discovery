# A script to ensure that technical replicates match for any given variant for any given DNA sample

library(magrittr)

# define new operator
`%notin%` = function(x,y) !(x %in% y)

# read in the data and reformat
square_vcf = readr::read_tsv(snakemake@input[['combined_vcfs']], comment='##')
square_vcf = square_vcf %>% dplyr::rename('CHROM'=`#CHROM`)

# remove the column in the table which corresponds to the normal sample
patient_id = snakemake@wildcards$patient_id %>% as.integer()

# remove the germline sample
germline_metadata = readr::read_tsv(snakemake@input[['germline_metadata']]) %>% dplyr::filter(fk_britroc_number==patient_id)
germline_sample = germline_metadata %>% dplyr::pull(fk_sample) %>% unique()
square_vcf = square_vcf %>% dplyr::select(-all_of(germline_sample))

# define samples and sample types
somatic_metadata = readr::read_tsv(snakemake@input[['matched_somatic_metadata']]) %>% dplyr::filter(fk_britroc_number==patient_id)
archival_samples = somatic_metadata %>% dplyr::filter(type=='archival') %>% dplyr::pull(fk_sample) %>% unique()
relapse_samples = somatic_metadata %>% dplyr::filter(type=='relapse') %>% dplyr::pull(fk_sample) %>% unique()
samples = c(archival_samples, relapse_samples)

# filter by quality score
square_vcf = square_vcf %>% dplyr::filter(QUAL>snakemake@params[['variant_quality_score_threshold']])

# write empty data frame to file if data frame is empty at this point and exit script
if (square_vcf %>% dim() %>% .[1] == 0) {
	readr::write_tsv(square_vcf %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['tumour_samples_union']])
	readr::write_tsv(square_vcf %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['archival_samples']])
	readr::write_tsv(sample_genotype_table %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['relapse_samples']])

	readr::write_tsv(square_vcf, path=snakemake@output[['library_MAFs']])
	readr::write_tsv(square_vcf, path=snakemake@output[['library_depths']])

	readr::write_tsv(square_vcf %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['sample_genotypes']])

	quit()
} else {
}

# TODO: formalise the process of mapping libraries to samples
identify_variants_with_tech_rep_mismatch_in_joined_vcf_table = function(square_vcf) {
	# square_vcf: A joined vcf table where each column represents one library
	# output - sample_genotype_table: For each variant in square_vcf, a genotype for each sample where multiple libraries map onto single samples

	# create space for creating a table of genotypes
	genotype_table = square_vcf

	# this assumed library names use sample names as a prefix (i.e. JBLAB-001 and JBLAB-001_d)
	libraries = genotype_table %>% dplyr::select(dplyr::starts_with(samples)) %>% colnames()

	# extract the predicted genotype for each variant
	for (lib in libraries) {
		genotype_table[[lib]] = genotype_table[[lib]] %>% stringr::str_extract(pattern='[01\\|]+')
	}

	# go from a library genotype table to a sample genotype table - first by removing all library fields
	sample_genotype_table = genotype_table %>% dplyr::select(-dplyr::any_of(libraries)) 

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

	for (lib in libraries) {
		square_vcf_ft[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[17]]

		genotype_table[[lib]] = ifelse(square_vcf_ft[[lib]] != 'PASS', '0|0', genotype_table[[lib]])
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
	sample_genotype_table = sample_genotype_table %>% dplyr::filter(dplyr::if_any(dplyr::all_of(samples), `%notin%`, c(NA,"0|0")))

	sample_genotype_table = sample_genotype_table %>% dplyr::inner_join(square_vcf %>% dplyr::select(CHROM,POS,REF,ALT), by=c('CHROM','POS','REF','ALT'))

	return(sample_genotype_table)
}

sample_genotype_table = identify_variants_with_tech_rep_mismatch_in_joined_vcf_table(square_vcf) 

# join the results of this test/filter into the main table
square_vcf = square_vcf %>% dplyr::inner_join(sample_genotype_table %>% dplyr::select(CHROM,POS,REF,ALT), by=c('CHROM','POS','REF','ALT'))

print(square_vcf, width=Inf)
quit()

#print(sample_genotype_table, n=5)

#print(square_vcf %>% unique(), n=5)

square_vcf = square_vcf %>% dplyr::group_by(CHROM,POS,REF,ALT) %>% dplyr::slice(n=1) %>% dplyr::ungroup()

#print(square_vcf %>% unique(), n=5)

sample_genotype_table = sample_genotype_table %>% dplyr::group_by(CHROM,POS,REF,ALT) %>% dplyr::slice(n=1) %>% dplyr::ungroup()

#print(sample_genotype_table %>% unique())

#print(square_vcf, width=Inf)

# get MAFs for each sample

#square_vcf$IM_403 %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[3]] #%>% 
	#stringr::str_split(pattern=',') %>% data.table::transpose() %>% .[[2]] %>% stringr::str_replace(pattern='^.$', replacement='0.000')

square_vcf_MAF = square_vcf
square_vcf_depth = square_vcf

#print(libraries)
#print(square_vcf_MAF)

#print('foo')

if (square_vcf_MAF %>% dim() %>% .[1] == 0) {
	readr::write_tsv(square_vcf_MAF %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['tumour_samples_union']])
	readr::write_tsv(square_vcf_MAF %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['archival_samples']])
	readr::write_tsv(square_vcf_MAF %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['relapse_samples']])

	readr::write_tsv(square_vcf_MAF, path=snakemake@output[['library_MAFs']])
	readr::write_tsv(square_vcf_depth, path=snakemake@output[['library_depths']])

	readr::write_tsv(sample_genotype_table, path=snakemake@output[['sample_genotypes']])

	quit()
} else {
	}


for (lib in libraries) {
	square_vcf_MAF[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[12]] %>%
		stringr::str_split(pattern=',') %>% data.table::transpose() %>% .[[2]] %>% stringr::str_replace(pattern='^.$', replacement='0.000') %>%
		as.numeric()
	square_vcf_depth[[lib]] = square_vcf[[lib]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[3]] %>%
		as.integer()
}

#print(square_vcf_MAF)
#print(square_vcf_depth)

sample_weighted_MAF = square_vcf %>% dplyr::select(-all_of(libraries))

for (sample in samples) {
	print(sample)
	relevant_depths = square_vcf_depth %>% dplyr::select(starts_with(sample))
	relevant_MAFs = square_vcf_MAF %>% dplyr::select(starts_with(sample))

	#relevant_libraries$sample_genotype = purrr::map(.x=1:dim(relevant_libraries)[1], .f=check_genotypes, relevant_libraries) %>% unlist()

	weighted_MAF_numerator = (relevant_depths[,1] * relevant_MAFs[,1]) + (relevant_depths[,2] * relevant_MAFs[,2])
	weighted_MAF_denominator = relevant_depths[,1] + relevant_depths[,2]
	weighted_MAF_numerator = weighted_MAF_numerator %>% dplyr::pull(sample)
	weighted_MAF_denominator = weighted_MAF_denominator %>% dplyr::pull(sample)

	maf_diff = abs(relevant_MAFs[,1] - relevant_MAFs[,2]) %>% dplyr::pull(sample)

	print(maf_diff)

	sample_weighted_MAF[[sample]] = weighted_MAF_numerator / weighted_MAF_denominator
	#sample_weighted_MAF = sample_weighted_MAF %>% dplyr::mutate(maf_diff = maf_diff)

	# ensure weighted mean MAF is above 5%
	
	#print('foo')

	sample_weighted_MAF = sample_weighted_MAF %>% tidyr::unite(col=substitution,REF,ALT, sep='>', remove=FALSE)

	#print(sample_weighted_MAF, width=Inf)
	#print(sample_genotype_table, width=Inf)

        #print(sample_genotype_table[[sample]])
	#print(sample_weighted_MAF[[sample]])

	new_sample_genotype_column = ifelse(
		sample_weighted_MAF[[sample]] < 0.23 & sample_weighted_MAF[['substitution']] == 'C>T',
		'0|0', 
		sample_genotype_table[[sample]]
		)

	#sample_genotype_table[[sample]] = new_sample_genotype_column

	new_sample_genotype_column = ifelse(
		sample_weighted_MAF[[sample]] < 0.23 & sample_weighted_MAF[['substitution']] == 'G>A',
		'0|0', 
		new_sample_genotype_column
		)

	#print(new_sample_genotype_column)

	new_sample_genotype_column = ifelse(
		maf_diff >= 0.30 & sample_weighted_MAF[['substitution']] == 'G>A',
		'0|0', 
		new_sample_genotype_column
		)

	#print(new_sample_genotype_column)

	new_sample_genotype_column = ifelse(
		maf_diff >= 0.30 & sample_weighted_MAF[['substitution']] == 'C>T',
		'0|0', 
		new_sample_genotype_column
		)

	#print(new_sample_genotype_column)
	#print(sample_genotype_table)
	#print(sample_weighted_MAF, width=Inf)
	#print(relevant_MAFs)

	#print(new_sample_genotype_column)

	#print('foo')
	
	sample_genotype_table[[sample]] = new_sample_genotype_column
	print(sample_genotype_table[[sample]])
}

sample_genotype_table = sample_genotype_table %>% dplyr::filter(dplyr::if_any(all_of(samples), `%notin%`, c(NA,"0|0")))

# filter original library genotype table
genotype_table = genotype_table %>% 
	dplyr::inner_join(sample_genotype_table %>% dplyr::select(-starts_with('IM'), -starts_with('JBLAB')), by=c('CHROM','POS','REF','ALT'))

#print(genotype_table %>% dplyr::select(starts_with('IM'), starts_with('JBLAB')), n=Inf)
#print(square_vcf_MAF %>% dplyr::select(starts_with('IM'), starts_with('JBLAB')), n=Inf)
#print(square_vcf_depth %>% dplyr::select(starts_with('IM'), starts_with('JBLAB')), n=Inf)
#print(sample_weighted_MAF %>% dplyr::select(starts_with('IM'), starts_with('JBLAB')), width=Inf, n=Inf)
#print(sample_genotype_table %>% dplyr::select(-FORMAT), width=Inf, n=Inf)

#print(sample_genotype_table, width=Inf)

# subset for variants which appear in at least one archival sample
sample_genotype_table_archival = sample_genotype_table %>% dplyr::filter(dplyr::if_any(all_of(archival_samples), `%notin%`, c(NA,"0|0")))

#print(sample_genotype_table_archival)

# subset for variants which appear in at least one relapse sample
sample_genotype_table_relapse = sample_genotype_table %>% dplyr::filter(dplyr::if_any(all_of(relapse_samples), `%notin%`, c(NA,"0|0")))

#print(square_vcf_MAF)
#print(square_vcf_depth)

#print(sample_genotype_table_relapse)

readr::write_tsv(square_vcf_MAF, path=snakemake@output[['library_MAFs']])
readr::write_tsv(square_vcf_depth, path=snakemake@output[['library_depths']])

readr::write_tsv(sample_genotype_table, path=snakemake@output[['sample_genotypes']])

readr::write_tsv(sample_genotype_table %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['tumour_samples_union']], append=FALSE)
readr::write_tsv(sample_genotype_table_archival %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['archival_samples']], append=FALSE)
readr::write_tsv(sample_genotype_table_relapse %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['relapse_samples']], append=FALSE)
