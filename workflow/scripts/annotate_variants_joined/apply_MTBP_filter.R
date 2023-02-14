library(magrittr)

input = readr::read_tsv(snakemake@input[[1]])
print(input)

MTBP_records = tibble::tribble(
	~CHROM, ~POS, ~REF, ~ALT,
	'chr17',	41245457,	'GA',	'G',
	'chr17',	41245586,	'CT',	'C',
	'chr17',	41243941,	'G',	'A',
	'chr17',	41234592,	'G',	'A',
	'chr13',	32907420,	'GA',	'G',
	'chr13',	32954022,	'CA',	'C',
	'chr13',	32914577,	'G',	'T',
	'chr13',	32913836,	'CA',	'C',
	'chr13',	32953632,	'CA',	'C',
	'chr13',	32913558,	'CA',	'C',
	'chr13',	32911442,	'GA',	'G',
	'chr2',		215646084,	'CT',	'C',
	'chr17',	41246839,	'C',	'A',
	'chr17',	41244132,	'CTTCCCATAGGCTG',	'C',
	'chr13',	32900668,	'TCTAGGAGCTGAGGTGGATC',	'T',				 					 	 
)

print(MTBP_records)

filtered_input = dplyr::semi_join(
		input,
		MTBP_records,
		by=c('CHROM','POS','REF','ALT')
	)

print(filtered_input)

filtered_input %>% dplyr::select(CHROM,POS,REF,ALT) %>% unique() %>% print()

readr::write_tsv(filtered_input, snakemake@output[[1]], col_names=TRUE)
