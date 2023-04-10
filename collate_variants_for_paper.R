# collate different variants types for paper


# unmatched germline ----

germline_variants = readr::read_tsv('results/variant_analysis/germline/octopus_unmatched/panel_6_28/merged/collated/filtered_vep_calls_octopus_joined.tsv')
germline_variants_haplotypecaller = readr::read_tsv('results/variant_analysis/germline/ampliconseq_pipeline/panel_6_28/collated/final_germline_variants.tsv')

#TODO: add additional germline information

# algorithm: haplotypecaller (add file)

# unmatched and unpaired ----

unmatched_unpaired_variants = readr::read_tsv('results/variant_analysis/unmatched/collated/filtered_vep_calls_octopus_joined.tsv')

# unmatched and paired ----

# corresponds to figure 2

unmatched_paired_archival = readr::read_tsv('results/variant_analysis/unmatched/paired/submission_results/collated/filtered_archival_vep_calls_octopus_joined.targeted.tsv')
unmatched_paired_relapse = readr::read_tsv('results/variant_analysis/unmatched/paired/submission_results/collated/filtered_relapse_vep_calls_octopus_joined.targeted.tsv')

# matched and unpaired ----

matched_and_unpaired_variants = readr::read_tsv('results/variant_analysis/matched/panel_6_28/collated/filtered_vep_calls_octopus_joined.tsv')

# matched and paired ----

matched_paired_archival = readr::read_tsv('results/variant_analysis/matched/panel_6_28/paired/collated/filtered_vep_calls_octopus_joined_archival.tsv')
matched_paired_relapse = readr::read_tsv('results/variant_analysis/matched/panel_6_28/paired/collated/filtered_vep_calls_octopus_joined_relapse.tsv')

# TP53 ----

TP53_archival = readr::read_tsv('results/variant_analysis/TP53/collated/filtered_archival_vep_calls_octopus.tsv')
TP53_relapse = readr::read_tsv('results/variant_analysis/TP53/collated/filtered_relapse_vep_calls_octopus.tsv')

unmatched_unpaired_variants = unmatched_unpaired_variants %>% 
	dplyr::select(-`#Uploaded_variation`,-Allele,-Feature,-Feature_type,
			-cDNA_position,-CDS_position,-Protein_position,-Amino_acids,-Codons,-Extra,-Location,
			-IMPACT,-Gene,-SYMBOL_SOURCE,-BIOTYPE,-SIFT,-PolyPhen,
			-CANONICAL,-ENSP,-UNIPARC,-EXON,
			-INTRON,-SWISSPROT,-TREMBL,-DOMAINS,
			-FLAGS)

matched_and_unpaired_variants = matched_and_unpaired_variants %>% 
	dplyr::select(-`#Uploaded_variation`,-Allele,-Feature,-Feature_type,
			-cDNA_position,-CDS_position,-Protein_position,-Amino_acids,-Codons,-Extra,-Location,
			-IMPACT,-Gene,-SYMBOL_SOURCE,-BIOTYPE,-SIFT,-PolyPhen,
			-CANONICAL,-ENSP,-UNIPARC,-EXON,
			-INTRON,-SWISSPROT,-TREMBL,-DOMAINS,
			-FLAGS)

germline_variants = germline_variants %>% 
	dplyr::select(-`#Uploaded_variation`,-Allele,-Feature,-Feature_type,
			-cDNA_position,-CDS_position,-Protein_position,-Amino_acids,-Codons,-Extra,-Location,
			-IMPACT,-Gene,-SYMBOL_SOURCE,-BIOTYPE,-SIFT,-PolyPhen,
			-CANONICAL,-ENSP,-UNIPARC,-EXON,
			-INTRON,-SWISSPROT,-TREMBL,-DOMAINS,
			-FLAGS,-ID,-sample_1,-sample_2,-QUAL,-INFO,-FORMAT,-FILTER)

germline_variants %>% colnames() %>% print()
germline_variants %>% dplyr::slice_head(n=1) %>% print(width=Inf)

germline_variants_haplotypecaller = germline_variants_haplotypecaller %>% dplyr::rename(
		'patient_id'='fk_britroc_number',
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
germline_variants_haplotypecaller = germline_variants_haplotypecaller %>% 
	dplyr::select(-fk_sample,-depth_1,-depth_2,-filter_1,-filter_2,-HGVS_united,-allele_fraction_1,-allele_fraction_2,-confidence,
	-united_effect,-fk_block_id,-fk_histological_id,-glasgow_comment)


germline_variants$sample_type = 'whole blood'
germline_variants$variant_caller = 'octopus'
germline_variants_haplotypecaller$analysis_type = 'unmatched_germline'
germline_variants_haplotypecaller$sample_type = 'whole blood'
germline_variants_haplotypecaller$variant_caller = 'HaplotypeCaller'

unmatched_unpaired_variants$sample_type = 'tumour'
unmatched_unpaired_variants$analysis_type = 'unmatched_and_unpaired'
unmatched_unpaired_variants$variant_caller = 'octopus'

matched_and_unpaired_variants$sample_type = 'tumour'
matched_and_unpaired_variants$analysis_type = 'matched_and_unpaired'
matched_and_unpaired_variants$variant_caller = 'octopus'

germline_variants$analysis_type = 'unmatched_germline'
#germline_variants$analysis_type = germline_variants$analysis_type %>% dplyr::recode(
#	'germline'='unmatched_germline'	
#)

#unmatched_unpaired_variants$fk_sample = NA
#unmatched_unpaired_variants$allele_fraction_1 = NA
#unmatched_unpaired_variants$allele_fraction_2 = NA

#setdiff(germline_variants %>% colnames(), unmatched_unpaired_variants %>% colnames())
#setdiff(unmatched_unpaired_variants %>% colnames(), germline_variants %>% colnames())

germline_variants$panel_id = 'panel_6_28_intersection'
germline_variants_haplotypecaller$panel_id = 'panel_6_28_intersection'
unmatched_unpaired_variants$panel_id = 'panel_28'
matched_and_unpaired_variants$panel_id = 'panel_6_28_intersection'

germline_variants %>% colnames() %>% print()
unmatched_unpaired_variants %>% colnames() %>% print()

germline_variants_clinical_records = tibble::tribble(
	~patient_id, ~Consequence, ~Existing_variation, ~CHROM, ~POS, ~REF, ~ALT, ~SYMBOL, ~CLIN_SIG, ~gnomAD_AF, ~AF, ~HGVSc, ~HGVSp, ~analysis_type, ~sample_type, ~variant_caller, ~panel_id,
	46, NA, NA, 'chr17', NA, NA, NA, 'BRCA1', NA, NA, NA, NA, NA, 'germline', 'whole blood', 'clinical testing', NA, 
	127, 'frameshift', NA, 'chr17', NA, NA, NA, 'BRCA1', NA, NA, NA, NA, 'p.Ala145_Asp1662fsX14', 'germline', 'whole blood', 'clinical testing', NA, 
	225, NA, NA, 'chr17', NA, NA, NA, 'BRCA1', NA, NA, NA, NA, NA, 'germline', 'whole blood', 'clinical testing', NA, 
	247, NA, NA, 'chr17', NA, NA, NA, 'BRCA1', NA, NA, NA, NA, NA, 'germline', 'whole blood', 'clinical testing', NA, 
	259, NA, NA, 'chr17', NA, NA, NA, 'BRCA1', NA, NA, NA, NA, NA, 'germline', 'whole blood', 'clinical testing', NA, 
	275, NA, NA, 'chr17', NA, NA, NA, 'BRCA1', NA, NA, NA, NA, NA, 'germline', 'whole blood', 'clinical testing', NA
)

print(germline_variants_clinical_records, width=Inf)

germline_variants_haplotypecaller %>% colnames() %>% print()

all_germline_variants  = rbind(
	germline_variants,
	germline_variants_clinical_records,
	germline_variants_haplotypecaller
	)


all_germline_variants = all_germline_variants %>% dplyr::filter(POS != '32954023')
all_germline_variants = all_germline_variants %>% dplyr::filter(!(patient_id == 145  & POS == '59885982'))
all_germline_variants = all_germline_variants %>% dplyr::filter(!(patient_id == 253  & POS == '32913551'))
all_germline_variants = all_germline_variants %>% dplyr::filter(!(patient_id == 274  & POS == '32913557'))

all_variants = rbind(
	all_germline_variants,
	unmatched_unpaired_variants,
	matched_and_unpaired_variants
	)
print('foo')

unmatched_paired_archival$sample_type = 'diagnosis tumour'
unmatched_paired_relapse$sample_type = 'relapse tumour'

matched_paired_archival$sample_type = 'diagnosis tumour'
matched_paired_relapse$sample_type = 'relapse tumour'

unmatched_paired_archival$variant_caller = 'octopus'
unmatched_paired_relapse$variant_caller = 'octopus'
matched_paired_archival$variant_caller = 'octopus'
matched_paired_relapse$variant_caller = 'octopus'

unmatched_paired_variants = rbind(unmatched_paired_archival, unmatched_paired_relapse)
matched_paired_variants = rbind(matched_paired_archival, matched_paired_relapse)
print('foo2')

unmatched_paired_variants$analysis_type = 'unmatched_paired'
matched_paired_variants$analysis_type = 'matched_paired'

unmatched_paired_variants$panel_id = 'panel_28'
matched_paired_variants$panel_id = 'panel_6_28_intersection'

matched_paired_variants %>% colnames() %>% print()
unmatched_paired_variants %>% colnames() %>% print()

paired_variants = rbind(matched_paired_variants, unmatched_paired_variants)

print('goo')

paired_variants = paired_variants %>% 
	dplyr::select(-`#Uploaded_variation`,-Allele,-Feature,-Feature_type,
			-cDNA_position,-CDS_position,-Protein_position,-Amino_acids,-Codons,-Extra,-Location,
			-IMPACT,-Gene,-SYMBOL_SOURCE,-BIOTYPE,-SIFT,-PolyPhen,
			-CANONICAL,-ENSP,-UNIPARC,-EXON,
			-INTRON,-SWISSPROT,-TREMBL,-DOMAINS,
			-FLAGS)

#unmatched_paired_variants$allele_fraction_1 = NA
#unmatched_paired_variants$allele_fraction_2 = NA

all_variants = 
	rbind(
		all_variants,
		paired_variants
	)

TP53_archival$sample_type = 'diagnosis tumour'
TP53_relapse$sample_type = 'relapse tumour'
TP53_archival$variant_caller = 'octopus'
TP53_relapse$variant_caller = 'octopus'

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
TP53_variants = TP53_variants %>% dplyr::select(-type, -`#Uploaded_variation`,-'classification',-sample_id)
TP53_variants = TP53_variants %>% unique()

TP53_variants$analysis_type = 'TP53'

TP53_variants = TP53_variants %>% dplyr::rename(
		'patient_id'='fk_britroc_number'
	)

TP53_variants = TP53_variants %>% tidyr::separate(Location, c('CHROM','POS'), sep=':', remove=TRUE)

#TP53_variants$allele_fraction_1 = NA
#TP53_variants$allele_fraction_2 = NA

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

print(all_variants %>% colnames())

readr::write_tsv(all_variants, 'short_variants_SI.tsv')
