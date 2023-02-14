library(magrittr)

input = readr::read_tsv(snakemake@input[[1]])
input = input %>% dplyr::select(-patient_id)

input = input %>% dplyr::mutate(
		ID = '.',
		QUAL = 5000,
		FILTER = 'PASS',
		INFO = 'foo'
	)

input = input %>% dplyr::select(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO)
input = input %>% dplyr::rename('#CHROM'='CHROM')
input = input %>% unique()

print(input)

readr::write_tsv(input, snakemake@output[[1]], col_names=TRUE)
