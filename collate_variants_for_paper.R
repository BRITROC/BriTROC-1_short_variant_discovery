# collate different variants types for paper

germline_variants = readr::read_tsv('results/variant_analysis/germline/panel_6_28/collated/final_germline_variants.tsv')

unmatched_unpaired_variants = readr::read_tsv('results/variant_analysis/unmatched/collated/filtered_vep_calls_octopus_joined.tsv')

unmatched_paired_archival = readr::read_tsv('results/variant_analysis/unmatched/collated/filtered_archival_vep_calls_octopus_joined.targeted.tsv')
unmatched_paired_relapse = readr::read_tsv('results/variant_analysis/unmatched/collated/filtered_relapse_vep_calls_octopus_joined.targeted.tsv')

TP53_archival = readr::read_tsv('results/variant_analysis/TP53/collated/filtered_archival_vep_calls_octopus.tsv')
TP53_relapse = readr::read_tsv('results/variant_analysis/TP53/collated/filtered_relapse_vep_calls_octopus.tsv')

#unmatched_paired_archival %>% colnames()
#unmatched_paired_relapse %>% colnames()
#TP53_archival %>% colnames()
#TP53_relapse %>% colnames()

germline_variants = germline_variants %>% dplyr::select(-confidence,-HGVSc,-fk_histological_id,-fk_block_id,-glasgow_comment)
unmatched_unpaired_variants = unmatched_unpaired_variants %>% 
	dplyr::select(-`#Uploaded_variation`,-Allele,-Feature,-Feature_type,
			-cDNA_position,-CDS_position,-Protein_position,-Amino_acids,-Codons,-Extra,-Location,
			-IMPACT,-Gene,-SYMBOL_SOURCE,-BIOTYPE,-SIFT,-PolyPhen,
			-CANONICAL,-ENSP,-UNIPARC,-EXON,
			-INTRON,-SWISSPROT,-TREMBL,-DOMAINS,
			-FLAGS)

germline_variants = germline_variants %>% dplyr::rename(
		'patient_id'='fk_britroc_number',
		'Gene'='gene_symbol',
		'Consequence'='variant_type',
		'HGVSc'='cdna_effect',
		'HGVSp'='protein_effect',
		'Existing_variation'='existing_variation',
		'CLIN_SIG'='clinical_signficance',
		'AF'='thousand_genomes_allele_frequency',
		'gnomAD_AF'='gnomad_allele_frequency',
		'analysis_type'='type',
		'REF'='ref',
		'ALT'='alt',
		'CHROM'='chromosome',
		'POS'='position',
		'SYMBOL'='gene_symbol'
	)

germline_variants$sample_type = 'non_tumour'
unmatched_unpaired_variants$sample_type = 'tumour'
unmatched_unpaired_variants$analysis_type = 'unmatched_and_unpaired'

germline_variants$analysis_type = germline_variants$analysis_type %>% dplyr::recode(
	'germline'='unmatched_germline'	
)

unmatched_unpaired_variants$fk_sample = NA
unmatched_unpaired_variants$allele_fraction_1 = NA
unmatched_unpaired_variants$allele_fraction_2 = NA

#setdiff(germline_variants %>% colnames(), unmatched_unpaired_variants %>% colnames())
#setdiff(unmatched_unpaired_variants %>% colnames(), germline_variants %>% colnames())

germline_variants$panel_id = 'panel_6_28_intersection'
unmatched_unpaired_variants$panel_id = 'panel_28'

all_variants = rbind(
	germline_variants,
	unmatched_unpaired_variants
	)

unmatched_paired_archival$sample_type = 'diagnosis tumour'
unmatched_paired_relapse$sample_type = 'relapse tumour'

unmatched_paired_variants = rbind(unmatched_paired_archival, unmatched_paired_relapse)
unmatched_paired_variants = unmatched_paired_variants %>% 
	dplyr::select(-`#Uploaded_variation`,-Allele,-Feature,-Feature_type,
			-cDNA_position,-CDS_position,-Protein_position,-Amino_acids,-Codons,-Extra,-Location,
			-IMPACT,-Gene,-SYMBOL_SOURCE,-BIOTYPE,-SIFT,-PolyPhen,
			-CANONICAL,-ENSP,-UNIPARC,-EXON,
			-INTRON,-SWISSPROT,-TREMBL,-DOMAINS,
			-FLAGS)
unmatched_paired_variants$analysis_type = 'unmatched_paired'
unmatched_paired_variants$fk_sample = NA
unmatched_paired_variants$allele_fraction_1 = NA
unmatched_paired_variants$allele_fraction_2 = NA

unmatched_paired_variants$panel_id = 'panel_28'

all_variants = 
	rbind(
		all_variants,
		unmatched_paired_variants
	)


TP53_archival$sample_type = 'diagnosis tumour'
TP53_relapse$sample_type = 'relapse tumour'

TP53_variants = rbind(TP53_archival, TP53_relapse)
TP53_variants = TP53_variants %>% 
	dplyr::select(-Allele,-Feature,-Feature_type,
			-cDNA_position,-CDS_position,-Protein_position,-Amino_acids,-Codons,
			-IMPACT,-Gene,-SYMBOL_SOURCE,-BIOTYPE,-SIFT,-PolyPhen,
			-CANONICAL,-ENSP,-UNIPARC,-EXON,
			-INTRON,-SWISSPROT,-TREMBL,-DOMAINS,
			-FLAGS)

TP53_variants$HGVSp = TP53_variants$Extra %>%
        stringr::str_extract('HGVSp=[A-Za-z0-9:\\.]+;') %>%
        stringr::str_remove(';$') %>%
        stringr::str_remove('HGVSp=[A-Za-z0-9\\.]+:')
TP53_variants$HGVSc = TP53_variants$Extra %>%
        stringr::str_extract('HGVSc=[A-Za-z0-9:\\.]+;') %>%
        stringr::str_remove(';$') %>%
        stringr::str_remove('HGVSc=[A-Za-z0-9\\.]+:')

TP53_variants = TP53_variants %>% dplyr::select(-Extra)

TP53_variants$tmp = TP53_variants$`#Uploaded_variation` %>% stringr::str_extract('_[-ACTG]+/[-ACTG]+') %>% stringr::str_remove('_')
TP53_variants = TP53_variants %>% tidyr::separate(`tmp`, c('REF','ALT'), sep='/', remove=TRUE)

TP53_driver_classification = readr::read_tsv('results/variant_analysis/TP53/collated/TP53_variants_with_clonality_classifications.tsv')
TP53_driver_classification = TP53_driver_classification %>% dplyr::select(fk_britroc_number, `#Uploaded_variation`,classification)

print(TP53_variants)
print(TP53_driver_classification)

TP53_variants = 
	dplyr::left_join(
		TP53_variants,
		TP53_driver_classification,
		by=c('fk_britroc_number','#Uploaded_variation')
)

TP53_variants = TP53_variants %>% dplyr::filter(classification=='clonal')
TP53_variants = TP53_variants %>% dplyr::select(-type, -`#Uploaded_variation`,-'classification')

TP53_variants$analysis_type = 'TP53'

TP53_variants = TP53_variants %>% dplyr::rename(
		'fk_sample'='sample_id',
		'patient_id'='fk_britroc_number'
	)

TP53_variants = TP53_variants %>% tidyr::separate(Location, c('CHROM','POS'), sep=':', remove=TRUE)

TP53_variants$allele_fraction_1 = NA
TP53_variants$allele_fraction_2 = NA

all_variants %>% colnames()
TP53_variants %>% colnames()

#setdiff(all_variants %>% colnames(), TP53_variants %>% colnames())
#setdiff(TP53_variants %>% colnames(), all_variants %>% colnames())

TP53_variants$panel_id = 'panel_1_10_28_union'

all_variants = 
	rbind(
		all_variants,
		TP53_variants
	)

readr::write_tsv(all_variants, 'short_variants_SI.tsv')
