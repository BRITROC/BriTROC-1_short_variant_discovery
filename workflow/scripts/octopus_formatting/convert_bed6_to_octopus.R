# convert genomic co-ordinates of amplicon targets from bed to octopus format

library(magrittr)

bed_intervals = readr::read_tsv(snakemake@input[[1]], col_names=FALSE)
colnames(bed_intervals) = c('chrom','start','end','amplicon','score','strand')

# sort bed file
bed_intervals$chrom_int = bed_intervals$chrom %>% stringr::str_remove('chr') %>% as.integer()

bed_intervals = bed_intervals %>% dplyr::arrange(chrom_int, start, end)

# select relevant columns
bed_intervals = bed_intervals %>% dplyr::select(chrom,start,end)

bed_intervals = bed_intervals %>% tidyr::unite(start,end, col='X_tmp', sep='-', remove=TRUE)

octopus_intervals = bed_intervals %>% tidyr::unite(chrom,X_tmp, col='X_final', sep=':', remove=TRUE)

readr::write_tsv(octopus_intervals, snakemake@output[[1]], col_names=FALSE)
