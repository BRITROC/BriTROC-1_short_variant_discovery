library(magrittr)
library(tidyverse)

# read in data
vcfs = readr::read_tsv('list_of_filtered_vcfs.txt', col_names=c('file_names'))

process_vcf_file = function(vcf_file_name) {

vcf_data = readr::read_tsv(vcf_file_name, comment='##', col_types=readr::cols(ALT = 'c', REF='c') )

# add sample name
if (vcf_data %>% dim() %>% .[1] != 0) {
	vcf_data$sample_name = stringr::str_extract(vcf_file_name, pattern='^(JBLAB-|IM_)[0-9]+')
} else { # account for samples with no variant calls
	vcf_data$sample_name =character() # add empty character vector to empty data frame to prevent downstream errors
}

print(vcf_data)

# select basic column as well as tumour FORMAT columns 
vcf_data = vcf_data %>% select(`#CHROM`,POS,REF,ALT,starts_with('IM'),starts_with('JBLAB'),sample_name)

# separate format columns for each tumour sample technical replicate
vcf_data = separate(vcf_data, col=6, into=c('GT_2','AD_2','AF_2','DP_2','F1R2_2','F2R1_2','SB_2'), sep=':') 
vcf_data = separate(vcf_data, col=5, into=c('GT_1','AD_1','AF_1','DP_1','F1R2_1','F2R1_1','SB_1'), sep=':') 

return(vcf_data)
}

x = map_dfr(.x=vcfs$file_names, .f=process_vcf_file) 

print(x)
print('foo')

# select relevant columns
x = x %>% select(sample_name,`#CHROM`,POS,REF,ALT,AD_1,AD_2) #AF_1,AF_2,DP_1,DP_2)

# separate AD field into different components
x = separate(x, col='AD_1', into=c('REF_AD_1','ALT_AD_1'), sep=',') 
x = separate(x, col='AD_2', into=c('REF_AD_2','ALT_AD_2'), sep=',') 

x$ALT_AD_1 = as.integer(x$ALT_AD_1)
x$REF_AD_1 = as.integer(x$REF_AD_1)
x$ALT_AD_2 = as.integer(x$ALT_AD_2)
x$REF_AD_2 = as.integer(x$REF_AD_2)

x$mutation_type = ifelse( length(x$REF) > 1 | length(x$ALT) > 1, 'indel', 'point mutation')

x$mutation_type = map2(.x=x$REF, .y=x$ALT, .f = function(x,y) if (nchar(x) > 1 | nchar(y) > 1) {return ('indel')} else {return ('point mutation')}     ) %>% unlist()

# calculate allele fractions
x$AF_1 = x$ALT_AD_1 / (x$ALT_AD_1 + x$REF_AD_1)
x$AF_2 = x$ALT_AD_2 / (x$ALT_AD_2 + x$REF_AD_2)

x$mean_AF = (x$AF_1 + x$AF_2) / 2

x$diff_AD = abs(x$AF_1 - x$AF_2)

# annotate variant with respect to identity in both technical replicates

x = x %>% filter(AF_1 >= 0.025 | AF_2 >= 0.025)

x$annotation = ifelse( x$ALT_AD_1 < 1 | x$ALT_AD_2 < 1, 'one rep only', 'both reps')

print(x %>% select(contains('REF'),contains('ALT'),annotation))

print(x$annotation %>% table())

x = x %>% filter(annotation=='both reps')

# filter for allele fraction
x = x %>% filter(AF_1 >= 0.025)
x = x %>% filter(AF_2 >= 0.025)

write_tsv(x, 'britroc1_mutect_high_confidence_calls.tsv')
