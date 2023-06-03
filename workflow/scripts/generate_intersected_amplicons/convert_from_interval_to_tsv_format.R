# A script to convert from interval to tsv format

library(magrittr)

amplicons = readr::read_tsv(
	snakemake@input[['amplicons']], 
	comment='@', 
	col_names=c('Chromosome', 'AmpliconStart','AmpliconEnd', 'Strand','ID')
	) %>% dplyr::select(-Strand)

targets = readr::read_tsv(
	snakemake@input[['targets']], 
	comment='@', 
	col_names=c('Chromosome', 'TargetStart','TargetEnd', 'Strand','ID')
	) %>% dplyr::select(-Strand)

merged_data = dplyr::inner_join(
		amplicons,
		targets,
		by=c('Chromosome','ID')
	)

merged_data = merged_data %>% dplyr::select(ID, Chromosome, AmpliconStart, AmpliconEnd, TargetStart, TargetEnd)
readr::write_tsv(merged_data, snakemake@output[[1]])
