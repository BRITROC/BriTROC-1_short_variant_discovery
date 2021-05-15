# A script to rejoin MAFs with ensembl annotated variants

library(magrittr)
library(ggplot2)
library(DBI)
library(RPostgres)
library(patchwork)

joined_calls = readr::read_tsv('foo.tsv')

#joined_calls %>% dplyr::group_by(sample_id,`#Uploaded_variation`) %>% dplyr::summarise(n=dplyr::n()) %>% dplyr::arrange(-n) %>% print()

#print(
#	joined_calls %>% dplyr::select(fk_britroc_number,sample_id,`#Uploaded_variation`,mean_AF) %>% dplyr::filter(is.na(mean_AF))
#)

# identify clonal and subclonal mutations
patient_level_calls = readr::read_tsv('foo2.tsv')

clonal_calls = dplyr::left_join(joined_calls,patient_level_calls, by=c('fk_britroc_number','#Uploaded_variation')) %>% dplyr::filter(classification=='clonal')

#print(clonal_calls %>% dplyr::select(fk_britroc_number,sample_id,`#Uploaded_variation`,mean_MAF), n=Inf)
# load in goranova_tp53_clonal_annotations

goranova_tp53_clonal = readr::read_tsv('resources/Goranova_tp53_clonal.tsv') %>% dplyr::select(SAMPLE_ID,TP53freq)

clonal_calls = dplyr::inner_join(clonal_calls,goranova_tp53_clonal, by=c('sample_id'='SAMPLE_ID'))

clonal_calls = clonal_calls %>% dplyr::filter(sample_id != 'IM_438')

#samples_with_no_call_by_tedi = clonal_calls %>% dplyr::filter(is.na(TP53freq)) %>% .$sample_id

#clonal_calls %>% dplyr::filter(is.na(TP53freq)) %>% dplyr::select(fk_britroc_number,sample_id,`#Uploaded_variation`,mean_MAF,TP53freq) %>% print()

#quit()

clonal_calls$TP53freq = ifelse(is.na(clonal_calls$TP53freq), 0, clonal_calls$TP53freq)

#clonal_calls %>% dplyr::filter(is.na(mean_MAF)) %>% print()

#clonal_calls %>% dplyr::filter(is.na(TP53freq)) %>% dplyr::select(fk_britroc_number,sample_id,mean_MAF,TP53freq) %>% print()

#quit()

#clonal_calls = clonal_calls %>% dplyr::filter(!is.na(mean_MAF))
#clonal_calls = clonal_calls %>% dplyr::filter(!is.na(TP53freq))

clonal_calls = clonal_calls %>% dplyr::mutate(MAF_diff = abs(mean_MAF - TP53freq)) %>% dplyr::arrange(-MAF_diff)

print(clonal_calls %>% dplyr::select('fk_britroc_number','sample_id','mean_MAF','TP53freq','MAF_diff'), n=Inf)

cor(clonal_calls$mean_MAF, clonal_calls$TP53freq)

#print(clonal_calls$`#Uploaded_variation` %>% unique())

clonal_calls$variant_type = ifelse(grepl('-',clonal_calls$`#Uploaded_variation`),'Short indel','SNV')
clonal_calls$variant_type = ifelse(grepl('CGGCTCATAGGGCACCACCACACTAT',clonal_calls$`#Uploaded_variation`),'Short indel', clonal_calls$variant_type)

clonal_calls %>% ggplot(aes(x=mean_MAF, y=TP53freq, colour=variant_type)) + geom_point() +
	geom_text(
		aes(label=ifelse(MAF_diff>0.15,as.character(sample_id),'')),
		hjust=0,
		vjust=0,
		#position=position_dodge(),
		colour='black',
		nudge_x = -0.035,
		nudge_y = 0.015,
		size=2
	) + xlab('TP53 MAF (T. Bradley)') + ylab('TP53 MAF (T. Goronova)') + ggtitle('BriTROC-1 TP53 driver mutations')

ggsave('TP53_clonal_calls_methods_comparison.png', device='png')

##### start proper analysis

clonal_calls = clonal_calls %>% dplyr::filter(MAF_diff < 0.15)

clonal_calls$Consequence = dplyr::recode(clonal_calls$Consequence,
	'inframe_deletion' = 'inframe indel',
	'inframe_insertion' = 'inframe indel',
	'stop_gained,frameshift_variant'='stop gained',
	'stop_gained,splice_region_variant'='stop gained',
	'stop_gained'='stop gained',
	'frameshift_variant'='frameshift variant',
	'missense_variant'='missense variant',
	'splice_acceptor_variant'='splice variant',
	'splice_donor_variant'='splice variant',
	'splice_region_variant,intron_variant'='splice variant',
	'splice_acceptor_variant,coding_sequence_variant,intron_variant'='splice variant'
)

#clonal_calls %>% dplyr::filter(grepl('splice_',Consequence)) %>% dplyr::select(`#Uploaded_variation`,Consequence) %>% print(width=Inf)

clonal_calls$Consequence = factor(clonal_calls$Consequence, levels=c('frameshift variant','stop gained','splice variant','inframe indel','missense variant'))

ggplot(clonal_calls, aes(x ='Consequence', fill=Consequence)) + geom_bar() + xlab('') + theme(
	legend.title=element_blank(),
	axis.text.x=element_blank()
	) + ggtitle('BriTROC-1: TP53 driver mutation type')

ggsave('tp53_mutation_type3.png', device='png')

####

colnames(clonal_calls) %>% print()

clonal_calls$DOMAINS = clonal_calls$DOMAINS %>% stringr::str_extract(pattern='Pfam:[A-Z0-9]+') %>% stringr::str_remove(pattern='Pfam:')
print(clonal_calls$DOMAINS %>% unique())

clonal_calls$DOMAINS = clonal_calls$DOMAINS %>% tidyr::replace_na('no domain identified')

clonal_calls$DOMAINS = dplyr::recode(clonal_calls$DOMAINS,
	'PF00870' = 'DNA-binding domain',
	'PF07710' = 'Tetramerisation domain',
)

p1 = ggplot(clonal_calls, aes(x =Consequence, fill=DOMAINS)) + geom_bar() + theme(
        legend.title=element_blank(),
        axis.title.x=element_blank()
        ) + ggtitle('BriTROC-1: p53 domains affected')
p2 = ggplot(clonal_calls, aes(x =Consequence, fill=DOMAINS)) + geom_bar(position='fill') + theme(                                                                                          legend.title=element_blank(),
        axis.title.x=element_blank()                                                                                                                                       ) + ggtitle('BriTROC-1: p53 domains affected') + ylab('Proportion')
p3 = p1 + p2

ggsave(p3, filename='tp53_domains4.png', device='png', scale=1.15, width=16)
### IHC

TMA_1_tp53 = readr::read_tsv('BriTROC-1_TMA_1_tp53.tsv')

TMA_1_tp53 = TMA_1_tp53 %>% dplyr::rename(sample_id='Sample_ID')
TMA_1_tp53 = TMA_1_tp53 %>% dplyr::select(sample_id,Mean_H_Score)

print(TMA_1_tp53)

TMA_2_tp53 = readr::read_tsv('BriTROC-1_TMA_2_tp53.tsv')

TMA_2_tp53 = TMA_2_tp53 %>% dplyr::rename(sample_id='IM_labID',Mean_H_Score='Mean')
TMA_2_tp53 = TMA_2_tp53 %>% dplyr::select(sample_id,Mean_H_Score)

TMA_tp53 = rbind(TMA_1_tp53,TMA_2_tp53)

print(TMA_tp53)

clonal_calls2 = dplyr::inner_join(clonal_calls,TMA_tp53, by='sample_id')

ggplot(clonal_calls2, aes(x=Consequence, y=Mean_H_Score)) +  geom_boxplot(aes(colour=Consequence)) + geom_point(aes(colour=Consequence), position = position_jitterdodge()) + theme(axis.title.x = element_blank(), legend.title=element_blank()) + ylab('p53 staining intensity ("mean H score")') +
	ggtitle('BriTROC-1: p53 expression by driver mutation type')

ggsave('tp53_IHC_boxplot.png', device='png')

######

signatures = readr::read_tsv('britroc_30kb_signature_data_sig_exposures.tsv')

signatures = t(signatures)
signatures = cbind(rownames(signatures), signatures)
colnames(signatures) = signatures[1,]
signatures = signatures[-1,]

signatures = signatures %>% tibble::as_tibble()

signatures = signatures %>% tidyr::pivot_longer(cols=c(s1,s2,s3,s4,s5,s6,s7), names_to='signature', values_to='signature_exposure')

print(signatures)

clonal_calls = dplyr::inner_join(clonal_calls, signatures, by=c('sample_id'='rowname'))

clonal_calls$signature = factor(clonal_calls$signature, levels=c('s1','s2','s3','s4','s5','s6','s7'))
clonal_calls$signature_exposure = as.numeric(clonal_calls$signature_exposure)

print(clonal_calls %>% dplyr::select(sample_id,Consequence,signature,signature_exposure))

ggplot(clonal_calls, aes(x=signature, y=signature_exposure)) + geom_boxplot(aes(colour=Consequence)) + geom_point(aes(colour=Consequence), position = position_jitterdodge())

ggsave('signatures_by_tp53_mutation22.png', width=32, device='png')

###

clonal_calls = clonal_calls %>% dplyr::filter(Consequence != 'splice variant')
clonal_calls$expression = clonal_calls$Consequence
clonal_calls$expression = clonal_calls$expression %>% dplyr::recode(
		'missense variant'='Likely GOF',
		'inframe indel'='Likely GOF',
		'stop gained'='Likely LOF',
		'frameshift variant'='Likely LOF'
)

#print(clonal_calls %>% 
#	dplyr::select(fk_britroc_number,signature,signature_exposure) %>% 
#	dplyr::filter(signature=='s1') %>%
#	dplyr::arrange(fk_britroc_number)
#)

one_sample_per_patient = clonal_calls %>% dplyr::group_by(fk_britroc_number) %>% dplyr::sample_n(1) %>% dplyr::ungroup() %>% .$sample_id

#clonal_calls3 = clonal_calls %>% dplyr::filter(sample_id %in% one_sample_per_patient)
#stages = c(1,2,3,4)
#sample_sizes = c(3,1,21,4)

#sample_from_group = function(FIGO_stage,sample_size) {
#	GOF = clonal_calls3 %>% dplyr::filter(expression=='Likely GOF') %>% dplyr::filter(stage==FIGO_stage) %>% dplyr::sample_n(sample_size)
#	LOF = clonal_calls3 %>% dplyr::filter(expression=='Likely LOF') %>% dplyr::filter(stage==FIGO_stage) %>% dplyr::sample_n(sample_size)
#	x = rbind(GOF,LOF)
#
#	return(x)
#	}

#clonal_calls3 = purrr::map2_dfr(.x=stages, .y=sample_sizes, .f=sample_from_group)

#ggplot(clonal_calls3,
#	 aes(x=signature, y=signature_exposure)
#	) + geom_boxplot(aes(colour=expression)) + geom_point(aes(colour=expression), position = position_jitterdodge())  + ggtitle('BriTROC-1: Copy number signature exposures vs. TP53 mutation type') + scale_colour_discrete(name = "TP53 mutation type") + 
#	guides(fill=guide_legend(title="TP53 mutation type")) 

#ggsave('signatures_by_tp53_expression.png', width=24, device='png')

ks.test(
	clonal_calls %>% dplyr::filter(expression=='Likely GOF') %>% dplyr::filter(signature=='s2') %>% .$signature_exposure,
	clonal_calls %>% dplyr::filter(expression=='Likely LOF') %>% dplyr::filter(signature=='s2') %>% .$signature_exposure,
)

ks.test(
	clonal_calls %>% dplyr::filter(expression=='Likely GOF') %>% dplyr::filter(signature=='s4') %>% .$signature_exposure,
	clonal_calls %>% dplyr::filter(expression=='Likely LOF') %>% dplyr::filter(signature=='s4') %>% .$signature_exposure,
)

britroc_con = dbConnect(RPostgres::Postgres(),
                 dbname='britroc1',
                 host='jblab-db.cri.camres.org', # use 'jblab-db.cruk.cam.ac.uk' instead for external IPs,
                 port = 5432,
                 user = Sys.getenv('jblab_db_username'),
                 password = Sys.getenv('jblab_db_password')
)

print(britroc_con)

patients = dbReadTable(britroc_con, 'patients') %>% dplyr::select(britroc_number,tumour_stage_at_diagnosis) %>%
	dplyr::rename('fk_britroc_number'='britroc_number','stage'='tumour_stage_at_diagnosis')
print(patients)

clonal_calls = dplyr::inner_join(clonal_calls, patients, by='fk_britroc_number')

#print(clonal_calls %>% dplyr::select(fk_britroc_number,stage,expression) %>% unique(), n=Inf)

clonal_calls$stage = factor(clonal_calls$stage, levels=c('1','2','3','4'))

p4 = ggplot(clonal_calls %>% dplyr::select(fk_britroc_number,stage,expression) %>% unique() %>% dplyr::filter(!is.na(stage)), 
	aes(x=expression, fill=stage)) + geom_bar() + xlab('TP53 mutation type') + ggtitle('BriTROC-1: original tumour stage versus TP53 mutation type')
p5 = ggplot(clonal_calls %>% dplyr::select(fk_britroc_number,stage,expression) %>% unique() %>% dplyr::filter(!is.na(stage)), 
	aes(x=expression, fill=stage)) + geom_bar(position='fill') + xlab('TP53 mutation type') + ylab('Proportion') + 
	ggtitle('BriTROC-1: original tumour stage versus TP53 mutation type')

p6 = p4 + p5

ggsave(p6, filename='tp53_tumour_stage4.png', device='png', width=16)

clonal_calls %>% dplyr::select(fk_britroc_number,stage,expression) %>% unique() %>% dplyr::filter(expression=='Likely GOF') %>% .$stage %>% table()
clonal_calls %>% dplyr::select(fk_britroc_number,stage,expression) %>% unique() %>% dplyr::filter(expression=='Likely LOF') %>% .$stage %>% table()

#######

clonal_calls3 = clonal_calls %>% dplyr::filter(sample_id %in% one_sample_per_patient)
stages = c(1,2,3,4)
sample_sizes = c(3,1,21,4)

sample_from_group = function(FIGO_stage,sample_size) {
	GOF_samples = clonal_calls3 %>% 
			dplyr::filter(expression=='Likely GOF') %>% 
			dplyr::filter(stage==FIGO_stage) %>% 
			dplyr::select(sample_id) %>% unique() %>% 
			dplyr::sample_n(sample_size) %>% .$sample_id
	
	LOF_samples = clonal_calls3 %>% 
			dplyr::filter(expression=='Likely LOF') %>%
			dplyr::filter(stage==FIGO_stage) %>% 
			dplyr::select(sample_id) %>% unique() %>% 
			dplyr::sample_n(sample_size) %>% .$sample_id

	GOF = clonal_calls3 %>% dplyr::filter(sample_id %in% GOF_samples)
	LOF = clonal_calls3 %>% dplyr::filter(sample_id %in% LOF_samples)
	x = rbind(GOF,LOF)

	return(x)
	}

clonal_calls3 = purrr::map2_dfr(.x=stages, .y=sample_sizes, .f=sample_from_group)

ggplot(clonal_calls3,
	 aes(x=signature, y=signature_exposure)
	) + geom_boxplot(aes(colour=expression)) + geom_point(aes(colour=expression), position = position_jitterdodge())  + ggtitle('BriTROC-1: Copy number signature exposures vs. TP53 mutation type') + scale_colour_discrete(name = "TP53 mutation type") + 
	guides(fill=guide_legend(title="TP53 mutation type")) 

ggsave('signatures_by_tp53_expression.png', width=24, device='png')

#####

germline_snvs = dbReadTable(britroc_con, 'germline_snvs') %>% tibble::as_tibble()
germline_indels = dbReadTable(britroc_con, 'germline_indels') %>% tibble::as_tibble()

germline_snvs = germline_snvs %>% dplyr::select(
			fk_sample,gene_symbol,allele_fraction_1,allele_fraction_2,variant_type,existing_variation,clinical_signficance,thousand_genomes_allele_frequency,gnomad_allele_frequency
		)
germline_indels =  germline_indels %>% dplyr::select(
                        fk_sample,gene_symbol,allele_fraction_1,allele_fraction_2,variant_type,existing_variation,clinical_signficance,thousand_genomes_allele_frequency,gnomad_allele_frequency                      )

germline_short_variants = rbind(germline_snvs,germline_indels)
germline_short_variants = germline_short_variants %>% dplyr::filter(gene_symbol %in% c('BRCA1','BRCA2')) %>% 
		dplyr::filter(thousand_genomes_allele_frequency < 0.001 | is.na(thousand_genomes_allele_frequency)) %>%
		dplyr::filter(gnomad_allele_frequency < 0.001 | is.na(gnomad_allele_frequency)) %>%
		dplyr::filter(!grepl('benign',clinical_signficance)) %>%
		dplyr::filter(! ( is.na(clinical_signficance) & variant_type %in% c('nonsynonymous','inframe_insertion','inframe_deletion'))) %>%
		dplyr::filter(! ( clinical_signficance=='uncertain_significance' & variant_type %in% c('nonsynonymous','inframe_insertion', 'inframe_deletion'))) %>% dplyr::filter(!variant_type %in% ('synonymous','intron','upstream_gene','downstream_gene','3_prime_UTR','5_prime_UTR'))

print(germline_short_variants, width=Inf, n=30)
