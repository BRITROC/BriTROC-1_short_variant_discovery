# convert genomic co-ordinates of amplicon targets from bed to octopus format

library(magrittr)

bed_intervals = readr::read_tsv(snakemake@input[[1]], col_names=FALSE)

# sort bed file
bed_intervals$chrom_int = bed_intervals$X1 %>% stringr::str_remove('chr') %>% as.integer()

bed_intervals = bed_intervals %>% dplyr::arrange(chrom_int, X7, X8)

# select relevant columns
bed_intervals = bed_intervals %>% dplyr::select(X1,X7,X8)

bed_intervals = bed_intervals %>% tidyr::unite(X7,X8, col='X_tmp', sep='-', remove=TRUE)

octopus_intervals = bed_intervals %>% tidyr::unite(X1,X_tmp, col='X_final', sep=':', remove=TRUE)

readr::write_tsv(octopus_intervals, snakemake@output[[1]], col_names=FALSE)
