# A script to ensure that technical replicates match for any given variant for any given DNA sample

library(magrittr)

# define new operator
`%notin%` = function(x,y) !(x %in% y)

# source functions
source('functions.R')

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

# this assumed library names use sample names as a prefix (i.e. JBLAB-001 and JBLAB-001_d)
libraries = square_vcf %>% dplyr::select(dplyr::starts_with(samples)) %>% colnames()

sample_genotype_table = identify_variants_with_tech_rep_mismatch_in_joined_vcf_table(square_vcf) 

# join the results of this test/filter into the main table
square_vcf = square_vcf %>% dplyr::inner_join(sample_genotype_table %>% dplyr::select(CHROM,POS,REF,ALT), by=c('CHROM','POS','REF','ALT'))

# remove duplicate recorded variants due to overlapping amplicons recording the same variant
square_vcf = square_vcf %>% dplyr::group_by(CHROM,POS,REF,ALT) %>% dplyr::slice(n=1) %>% dplyr::ungroup()
sample_genotype_table = sample_genotype_table %>% dplyr::group_by(CHROM,POS,REF,ALT) %>% dplyr::slice(n=1) %>% dplyr::ungroup()

# get MAFs for each sample
square_vcf_MAF = square_vcf
square_vcf_depth = square_vcf

# write straight to disk if variant table is empty at this point
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

sample_genotype_table = implement_substitution_type_specific_filters(square_vcf)

# subset for variants which appear in at least one archival sample
sample_genotype_table_archival = sample_genotype_table %>% dplyr::filter(dplyr::if_any(all_of(archival_samples), `%notin%`, c(NA,"0|0")))

# subset for variants which appear in at least one relapse sample
sample_genotype_table_relapse = sample_genotype_table %>% dplyr::filter(dplyr::if_any(all_of(relapse_samples), `%notin%`, c(NA,"0|0")))

# write data objects to disk
readr::write_tsv(square_vcf_MAF, path=snakemake@output[['library_MAFs']])
readr::write_tsv(square_vcf_depth, path=snakemake@output[['library_depths']])
readr::write_tsv(sample_genotype_table, path=snakemake@output[['sample_genotypes']])
readr::write_tsv(sample_genotype_table %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['tumour_samples_union']], append=FALSE)
readr::write_tsv(sample_genotype_table_archival %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['archival_samples']], append=FALSE)
readr::write_tsv(sample_genotype_table_relapse %>% dplyr::select(CHROM,POS,REF,ALT), path=snakemake@output[['relapse_samples']], append=FALSE)
