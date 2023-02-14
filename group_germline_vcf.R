# A realtively simple script which groups VCFs in small numbers of no greater than 20 in order to perform joint variant calling

list_of_vcfs = list.files(path='results/variant_analysis/germline', pattern='*.passed.vcf', full.names=TRUE)

# construct a tibble from the list of VCFs
# and assign each record to groups of no more than 20
vcf_table = tibble::tibble(vcf_name=list_of_vcfs, bam_name=list_of_vcfs, group=rep_len(1:20, length.out=length(list_of_vcfs)))
vcf_table$SLX_ID = stringr::str_extract(vcf_table$vcf_name, pattern='SLX-[0-9]+')

# correct bam path
vcf_table$bam_name = vcf_table$bam_name %>% stringr::str_replace('results/variant_analysis/germline/','../SLX/SLX_ID/bam/cleaned_bams/')
vcf_table$bam_name = vcf_table$bam_name %>% stringr::str_replace('SLX_ID',vcf_table$SLX_ID)
vcf_table$bam_name = vcf_table$bam_name %>% stringr::str_replace('passed.vcf','bam')

vcf_table$SLX_ID = NULL

print(vcf_table)

readr::write_tsv(vcf_table, file='results/variant_analysis/germline/vcf_table.tsv')
