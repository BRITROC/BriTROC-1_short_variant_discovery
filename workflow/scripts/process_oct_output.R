# A script for processing output from the octopus variant caller

library(tidyverse)

vcf = readr::read_tsv(snakemake@input[[1]], comment='##') #'results/tumour_sample_vcfs_octopus/JBLAB-19326.vcf', comment='##')  #JBLAB-4193_octopus.vcf', comment='##')

# filter for tp53 region
vcf = vcf %>% dplyr::filter(`#CHROM`=='chr17') %>% dplyr::filter(POS>7000000) %>% dplyr::filter(POS<8000000) 
# filter for PASS
#vcf = vcf %>% dplyr::filter(FILTER=='PASS') 
# filter for SOMATIC
vcf = vcf %>% dplyr::filter(grepl('SOMATIC',INFO))

#print(vcf$INFO)
#print(vcf$FORMAT)
#print(vcf$`SLX-14363.FLD0001`)
#print(vcf$`SLX-14363.FLD0025`)

print(vcf %>% colnames())

vcf = tidyr::separate(vcf, snakemake@wildcards[['sample']], into=c('GT_1','GQ_1','DP_1','MQ_1','PS_1','PQ_1','HSS_1','HPC_1','MAP_HF_1','HF_CR_1','FT_1'), sep=':', remove=TRUE)
vcf = tidyr::separate(vcf, paste(snakemake@wildcards[['sample']],'d',sep='_'), into=c('GT_2','GQ_2','DP_2','MQ_2','PS_2','PQ_2','HSS_2','HPC_2','MAP_HF_2','HF_CR_2','FT_2'), sep=':', remove=TRUE)
vcf = vcf %>% dplyr::select(-ID,-FILTER)

# filter so the genotypes of the two technical replicates match as well as the phase set
#vcf = vcf %>% dplyr::filter(GT_1==GT_2)
#vcf = vcf %>% dplyr::filter(PS_1==PS_2)

# filter for quality score
#vcf = vcf %>% dplyr::filter(QUAL<10000)

print(vcf %>% colnames)

print(vcf %>% select(INFO))
print(vcf %>% select(FORMAT))

print(vcf %>% select(POS,REF,ALT, QUAL, GT_1, GT_2, PS_1,MAP_HF_1, MAP_HF_2), width=Inf, n=Inf)  #MAP_HF_2))    #MAP_HF_2))
