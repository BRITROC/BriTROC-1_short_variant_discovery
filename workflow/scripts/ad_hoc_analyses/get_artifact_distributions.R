
library(magrittr)
library(ggplot2)
library(patchwork)

vcf = readr::read_tsv('results/tumour_sample_vcfs_octopus/39.filtered2.vcf', comment='##') #%>%
	#dplyr::select(JBLAB-4147   JBLAB-4147_d JBLAB-4147   JBLAB-4147_d JBLAB-4147   JBLAB-4147_d IM_17   IM_17_d JBLAB-4147      JBLAB-4147_d)

vcf$`JBLAB-4147_GT` = vcf$`JBLAB-4147` %>% stringr::str_extract('^[^:]+')
vcf$`JBLAB-4147_d_GT` = vcf$`JBLAB-4147_d` %>% stringr::str_extract('^[^:]+')

# get calls which only appear in one of the replicates

num_genuine_variants = vcf %>% dplyr::filter(`JBLAB-4147_GT` == `JBLAB-4147_d_GT`) %>% dplyr::filter(`JBLAB-4147_GT` != '0|0')

vcf = vcf %>% dplyr::filter(`JBLAB-4147_GT` != '0|0' | `JBLAB-4147_d_GT` != '0|0' )
vcf = vcf %>% dplyr::filter(`JBLAB-4147_GT` !=  `JBLAB-4147_d_GT` )

# determine the substitution type
vcf = vcf %>% tidyr::unite(substitution, REF, ALT, remove=TRUE, sep='>')

# get depth
vcf$depth_1 = vcf$`JBLAB-4147` %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[3]] %>%
                as.integer()
vcf$depth_2 = vcf$`JBLAB-4147_d` %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[3]] %>%
                as.integer()

vcf$FT_1 = vcf$`JBLAB-4147` %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[17]]
vcf$FT_2 = vcf$`JBLAB-4147_d` %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[17]]

vcf$MAF_1 = vcf$`JBLAB-4147` %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[12]] %>%
                stringr::str_split(pattern=',') %>% data.table::transpose() %>% .[[2]] %>% stringr::str_replace(pattern='^.$', replacement='0.000') %>%
                as.numeric()
vcf$MAF_2 = vcf$`JBLAB-4147_d` %>% stringr::str_split(pattern=':') %>% data.table::transpose() %>% .[[12]] %>%
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

lib_1 = vcf %>% dplyr::filter(`JBLAB-4147_GT` != '0|0') %>% dplyr::filter(depth_2 > 75) %>% dplyr::filter(FT_1 == 'PASS')
lib_2 = vcf %>% dplyr::filter(`JBLAB-4147_d_GT` != '0|0') %>% dplyr::filter(depth_1 > 75) %>% dplyr::filter(FT_2 == 'PASS')

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
p = p + plot_annotation(title='JBLAB-4147')

ggsave(plot=p, 'plots/JBLAB-4147_C_to_T.png', device='png', width=32)
#ggsave(plot=substitution_plot, 'plots/JBLAB-4147_substitution.png'  ,device='png')

#ggsave(plot=lib_2_plot, 'plots/JBLAB-4147_d_C_to_T.png', device='png')



