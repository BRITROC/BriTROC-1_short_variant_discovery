# A script to filter output variants based on curation by MTBP

library(magrittr)

# MTBP usage details
#	oncoKB access token used
#	vcf input
#	analysis run date: 05/02/2023 15:25
#	pipeline version used: 7.2.0

filtered_germline_variants = readr::read_tsv(snakemake@input[[1]])
print(filtered_germline_variants)

MTBP_variants_to_keep = tibble::tribble(
	~chromosome, ~position, ~ref, ~alt,
	'chr17',	41244865,	'GTT',	'G',
	'chr17',	41228595,	'AT',	'A',
	'chr17',	41243788,	'TAGAC',	'T',
	'chr17',	41243876,	'G',	'GGGAA',
	'chr17',	41243941,	'G',	'A',
	'chr17',	41244318,	'CCT',	'C',
	'chr17',	41244542,	'GT',	'G',
	'chr17',	41244907,	'C',	'A',
	'chr17',	41226499,	'C',	'T',
	'chr17',	41245072,	'TG',	'T',
	'chr17',	41245354,	'C',	'A',
	'chr17',	41245473,	'TG',	'T',
	'chr17',	41245792,	'GTTCAGCT',	'G',
	'chr17',	41246651,	'TAC',	'T',
	'chr17',	41267755,	'T',	'C',
	'chr17',	41276044,	'ACT',	'A',
	'chr17',	41226515,	'G',	'T',
	'chr17',	41223063,	'G',	'C',
	'chr17',	41223056,	'A',	'C',
	'chr17',	41215920,	'G',	'T',
	'chr17',	41209139,	'A',	'G',
	'chr17',	41226540,	'T',	'C',
	'chr17',	41256279,	'CT',	'C',
	'chr13',	32890621,	'GC',	'G',
	'chr13',	32914437,	'GT',	'G',
	'chr13',	32903604,	'CTG',	'C',
	'chr13',	32906804,	'C',	'CTTAG',
	'chr13',	32911595,	'G',	'T',
	'chr13',	32912895,	'CTGACA',	'C',
	'chr13',	32912964,	'TGAAA',	'T',
	'chr13',	32913199,	'CAG',	'C',
	'chr13',	32913558,	'C',	'CA',
	'chr13',	32913836,	'CAA',	'C',
	'chr13',	32914766,	'CTT',	'C',
	'chr13',	32929349,	'GA',	'G',
	'chr13',	32944692,	'C',	'T',
	'chr13',	32954022,	'CA',	'C',
	'chr13',	32913550,	'ACTTG',	'A',
	'chr13',	32913556,	'AGC',	'A',
	'chr16',	23625409,	'AT',	'A',
	'chr14',	68292235,	'C',	'T',
	'chr17',	56772543,	'C',	'T',
	'chr17',	33428320,	'C',	'T',
	'chr17',	59885981,	'CTGCT',	'C'		
)

final_germline_variants = dplyr::semi_join(
		filtered_germline_variants,
		MTBP_variants_to_keep,
		by=c('chromosome','position','ref','alt')		
	)

final_germline_variants %>% dplyr::select(chromosome,position,ref,alt) %>% unique() %>% print()

readr::write_tsv(final_germline_variants, file=snakemake@output[[1]])

quit()
