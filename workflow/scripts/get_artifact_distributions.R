
library(magrittr)
library(ggplot2)
library(patchwork)

vcf = readr::read_tsv(snakemake@input[[1]], comment='##') #%>%
	#dplyr::select(JBLAB-4147   JBLAB-4147_d JBLAB-4147   JBLAB-4147_d JBLAB-4147   JBLAB-4147_d IM_17   IM_17_d JBLAB-4147      JBLAB-4147_d)

sample_id = snakemake@wildcards$sample
tech_rep_2_id = paste(sample_id, 'd', sep='_')

print(sample_id)


vcf$`GT_tech_rep_1` = vcf[[sample_id]] %>% stringr::str_extract('^[^:]+')

vcf$`GT_tech_rep_2` = vcf[[tech_rep_2_id]] %>% stringr::str_extract('^[^:]+')


# get calls which only appear in one of the replicates

num_genuine_variants = vcf %>% dplyr::filter(`GT_tech_rep_1` == `GT_tech_rep_2`) %>% dplyr::filter(`GT_tech_rep_1` != '0|0')

vcf = vcf %>% dplyr::filter(`GT_tech_rep_1` != '0|0' | `GT_tech_rep_2` != '0|0' )
vcf = vcf %>% dplyr::filter(`GT_tech_rep_1` !=  `GT_tech_rep_2` )

# determine the substitution type
vcf = vcf %>% tidyr::unite(substitution, REF, ALT, remove=TRUE, sep='>')

# get depth
vcf$depth_1 = vcf[[sample_id]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[3]] %>%
                as.integer()
vcf$depth_2 = vcf[[tech_rep_2_id]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[3]] %>%
                as.integer()

vcf$FT_1 = vcf[[sample_id]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[17]]
vcf$FT_2 = vcf[[tech_rep_2_id]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[17]]

vcf$MAF_1 = vcf[[sample_id]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[12]] %>%
                stringr::str_split(pattern=',') %>% data.table::transpose() %>% .[[2]] %>% stringr::str_replace(pattern='^.$', replacement='0.000') %>%
                as.numeric()
vcf$MAF_2 = vcf[[tech_rep_2_id]] %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[12]] %>%
                stringr::str_split(pattern=',') %>% data.table::transpose() %>% .[[2]] %>% stringr::str_replace(pattern='^.$', replacement='0.000') %>%
                as.numeric()

print(vcf, width=Inf)

vcf$substitution = dplyr::recode(vcf$substitution,
                                               'G>T'='C>A',
                                               'G>C'='C>G',
                                               'G>A'='C>T',
                                               'A>T'='T>A',
                                               'A>G'='T>C',
                                               'A>C'='T>G'
                                               )

lib_1 = vcf %>% dplyr::filter(`GT_tech_rep_1` != '0|0') %>% dplyr::filter(depth_2 > 75) %>% dplyr::filter(FT_1 == 'PASS')
lib_2 = vcf %>% dplyr::filter(`GT_tech_rep_2` != '0|0') %>% dplyr::filter(depth_1 > 75) %>% dplyr::filter(FT_2 == 'PASS')

print(
	vcf %>% dplyr::arrange(-MAF_2) %>% dplyr::filter(FT_1 == 'PASS' & FT_2 == 'PASS') %>% dplyr::filter(substitution == 'C>T'), 
	width=Inf, n=Inf
)
lib_1$substitution %>% table()
lib_2$substitution %>% table()

print('num_genuine_variants')
print(num_genuine_variants %>% dim() %>% .[1])

lib_1$MAF_1 %>% summary()
lib_2$MAF_2 %>% summary()

lib_1$rep_num = 'tech_rep_1'
lib_2$rep_num = 'tech_rep_2'

lib_1$MAF = lib_1$MAF_1
lib_2$MAF = lib_2$MAF_2

lib = rbind(lib_1, lib_2)

print(lib)

lib_plot = ggplot(lib %>% dplyr::filter(substitution=='C>T'),  aes(x=MAF, fill=rep_num, colour=rep_num)  ) + geom_density(alpha=0.2) + geom_point(y=0) + xlim(c(0,1)) +
	ggtitle('C>T substitutions')
#lib_2_plot = ggplot(lib_2 %>% dplyr::filter(substitution=='C>T'),  aes(x=`MAF_2`)  ) + geom_density() + geom_point(y=0)

substitution_plot =  ggplot(lib, aes(x=substitution, fill=rep_num)) + geom_bar(position = position_dodge2(preserve = "single")) + ylim(c(0,120)) +
	ggtitle('Substitution type')

p = lib_plot | substitution_plot
p = p + plot_annotation(title=sample_id)

print(lib)

summary = tibble::tibble(sample_id=sample_id, num_artifacts=dim(lib) %>% .[1], MAF_upper_quartile=quantile(lib$MAF, c(0.95)))

ggsave(plot=p, snakemake@output[[1]], device='png', width=32)
readr::write_tsv(summary, snakemake@output[[2]])
#ggsave(plot=substitution_plot, 'plots/JBLAB-4147_substitution.png'  ,device='png')

#ggsave(plot=lib_2_plot, 'plots/JBLAB-4147_d_C_to_T.png', device='png')



