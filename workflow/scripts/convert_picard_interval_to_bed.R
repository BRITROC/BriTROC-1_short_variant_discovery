# A simple script to convert picard intervals formatted file to a specific 8 field BED format used by VarDict. This 8 field format contains genomic intervals for the amplicon in fields 3 and 4, and genomic intervals for the amplicon targeted region (i.e. amplicon without primers) in fields 7 and 8.
# The bed file produced is also used as a generic genomic interval resource for producing genomic interval data in octopus format

library(magrittr)

picard_interval_targets_file = readr::read_tsv(snakemake@input[['targets']], comment='@', col_names=c('chromosome','start','end','strand','amplicon'))
picard_interval_amplicons_file = readr::read_tsv(snakemake@input[['amplicons']], comment='@', col_names=c('chromosome','start','end','strand','amplicon'))

print(picard_interval_targets_file)
print(picard_interval_amplicons_file)

picard_interval_joined = dplyr::inner_join(
	picard_interval_amplicons_file, 
	picard_interval_targets_file, 
	by=c('chromosome','strand','amplicon'), suffix=c('_amplicons','_targets')
) %>% dplyr::mutate(
	score='.'
)

bed_file = picard_interval_joined %>% dplyr::select(chromosome, start_amplicons, end_amplicons, amplicon, score, strand, start_targets, end_targets)

readr::write_tsv(bed_file, snakemake@output[[1]], col_names=FALSE)
