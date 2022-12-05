# process a VCF of germline variants from a matched analysis in order to remove those variants which did not pass QC

library(magrittr)

all_germline_variants = readr::read_tsv(snakemake@input[[1]], comment='##')

# retrieve germline sample ID for this patient ----
germline_metadata = readr::read_tsv('config/germline_metadata_panel_matched_and_unpaired.tsv')
germline_sample_name = germline_metadata %>% dplyr::filter(fk_britroc_number==snakemake@wildcards['patient_id']) %>% .$fk_sample %>% unique()

# select relevant columns ----
#all_germline_variants = all_germline_variants %>% dplyr::select(`#CHROM`,POS,ALT,REF,ALT,QUAL,FORMAT,germline_sample_name)
#print(all_germline_variants, width=Inf)

# select only those germline variants which passed filters ----
all_germline_variants = all_germline_variants %>% dplyr::filter(grepl(':PASS$',x=!!as.name(germline_sample_name) ))
print(all_germline_variants, width=Inf)

readr::write_tsv(all_germline_variants, snakemake@output[[1]], col_names=TRUE)
