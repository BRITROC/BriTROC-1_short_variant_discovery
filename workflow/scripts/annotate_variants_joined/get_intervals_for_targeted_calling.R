# A script which takes filtered calls from non-targeted octopus variant calls and produces a set of intervals for targeted calling
# The idea is to test the set of intervals generated in one set of samples in another set of samples

library(magrittr)

filter3_calls = readr::read_tsv(snakemake@input[[1]])

print(filter3_calls)

# generate a column corresponding to the correct octopus interval for each record

# The interval pattern will be different depending on whether the record is an SNV or an indel
# identify the maximum allele length
output = filter3_calls %>% dplyr::mutate(
	MAX_ALLELE_LENGTH=dplyr::if_else(nchar(REF)>nchar(ALT), nchar(REF), nchar(ALT))
)

# construct interval
output = output %>% dplyr::mutate(
	INTERVAL=dplyr::if_else(MAX_ALLELE_LENGTH==1, as.character(POS-1), paste(POS-MAX_ALLELE_LENGTH,POS+MAX_ALLELE_LENGTH,sep='-'))
)

# prepend chromosome information
output = output %>% dplyr::mutate(
	INTERVAL=paste(CHROM,INTERVAL,sep=':')
)

# select only the interval column
output = output %>% dplyr::select(INTERVAL)

# write to file
readr::write_tsv(output, snakemake@output[[1]], col_names=FALSE)
