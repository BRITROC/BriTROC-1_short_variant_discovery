# A script to compare unmatched germline results for ampliconseq and octopus

library(DBI)
library(RPostgres)

# read-in data ----

source('functions.R')
source('~/.Renviron')

ampliconseq_germline = readr::read_tsv('results/variant_analysis/germline/ampliconseq_pipeline/panel_6_28/collated/final_germline_variants.tsv')
octopus_germline = readr::read_tsv('results/variant_analysis/germline/octopus_unmatched/panel_6_28/merged/collated/octopus_unmatched_germline_collated_final.vcf')

# remove amplicon duplicates ----

ampliconseq_germline = 
	ampliconseq_germline %>% 
	dplyr::group_by(fk_britroc_number,chromosome,position,ref,alt) %>%
	dplyr::slice_head(n=1) %>%
	dplyr::ungroup()

octopus_germline = 
	octopus_germline %>% 
	dplyr::group_by(patient_id,`#CHROM`,POS,REF,ALT) %>%
	dplyr::slice_head(n=1) %>%
	dplyr::ungroup()

# standardise the reporting of indel positions

ampliconseq_germline$position = 
	dplyr::if_else(
		nchar(ampliconseq_germline$ref) > 1 | nchar(ampliconseq_germline$alt) > 1,
		ampliconseq_germline$position + 1,
		ampliconseq_germline$position
	)

# filter

#ampliconseq_germline = ampliconseq_germline %>% dplyr::filter(allele_fraction_1>=0.20)
ampliconseq_germline = ampliconseq_germline %>% dplyr::filter(!(confidence=='medium' & !is.na(allele_fraction_2)))

ampliconseq_germline$HGVS_united %>% table() %>% sort()

# find overlap between sets

ampliconseq_exclusive = 
	dplyr::anti_join(
		ampliconseq_germline,
		octopus_germline,
		by=c('fk_britroc_number'='patient_id', 'chromosome'='#CHROM', 'position'='POS', 'ref'='REF', 'alt'='ALT')
	)

print(ampliconseq_exclusive)

octopus_exclusive = 
	dplyr::anti_join(
		octopus_germline,
		ampliconseq_germline,
		by=c('patient_id'='fk_britroc_number', '#CHROM'='chromosome', 'POS'='position', 'REF'='ref', 'ALT'='alt')
	)

print(octopus_exclusive)

ampliconseq_germline %>% dplyr::filter(fk_britroc_number==120)
octopus_germline %>% dplyr::filter(patient_id==120) %>% print(width=Inf)

ampliconseq_germline %>% dplyr::arrange(allele_fraction_1) %>% print(width=Inf)

#quit()

# generate a final table to report results ----

octopus_germline = octopus_germline %>%
	dplyr::rename(
		fk_britroc_number = patient_id,
		chromosome = `#CHROM`,
		position = POS,
		ref = REF,
		alt = ALT	
	)

ampliconseq_germline_reduced = ampliconseq_germline %>% 
	dplyr::select(fk_britroc_number,gene_symbol,chromosome,position,ref,alt,variant_type,HGVS_united) %>%
	dplyr::mutate(called_by=TRUE)
octopus_germline_reduced = octopus_germline %>% 
	dplyr::select(fk_britroc_number, chromosome, position, ref, alt) %>%
	dplyr::mutate(called_by=TRUE)

print(ampliconseq_germline_reduced)
print(octopus_germline_reduced)

combined_results = dplyr::full_join(
		ampliconseq_germline_reduced,
		octopus_germline_reduced,
		by=c('fk_britroc_number','chromosome','position','ref','alt'),
		suffix=c('_HaplotypeCaller','_octopus')
	)

combined_results$called_by_HaplotypeCaller = combined_results$called_by_HaplotypeCaller %>% tidyr::replace_na(FALSE)
combined_results$called_by_octopus = combined_results$called_by_octopus %>% tidyr::replace_na(FALSE)

combined_results$germline_sequencing = TRUE

combined_results$final_call_set = dplyr::if_else(
	combined_results$called_by_octopus == TRUE, # this is because octopus correctly calls a clinical record which haplotypecaller missed
	TRUE,
	FALSE
	)

combined_results$final_call_set = dplyr::if_else(
	combined_results$fk_britroc_number == 191, # A HaplotypeCaller exclusive call corroborated by clinical records - not called by octopus due to strand bias
	TRUE,
	combined_results$final_call_set
	)

combined_results$comment = dplyr::if_else(
	combined_results$fk_britroc_number == 120,
	'Variant corroborated by clinical reports so variant preserved in final call set',
	as.character(NA)
	)

combined_results$comment = dplyr::if_else(
	combined_results$fk_britroc_number == 91,
	'Patient BRCA1 status corroborated by clinical reports so variant preserved in final call set',
	as.character(NA)
	)

combined_results$final_call_set = dplyr::if_else(
	combined_results$HGVS_united == 'BRCA2:p.T3033Lfs*29',
	FALSE,
	combined_results$final_call_set,
	combined_results$final_call_set
	)

combined_results$comment = dplyr::if_else(
	combined_results$HGVS_united == 'BRCA2:p.T3033Lfs*29',
	'Low MAF; HGVS for some cases fails QC filter; does not pass QC filter in tumour samples; low plausibility double BRCA1/BRCA2 mutation for case 120',
	combined_results$comment,
	combined_results$comment
	)

combined_results$HGVS_united = combined_results$HGVS_united %>% stringr::str_remove('[A-Z0-9]+:')
combined_results$variant_type = combined_results$variant_type %>% dplyr::recode('nonsynonymous'='missense')

combined_results$gene_symbol  = combined_results$gene_symbol %>% tidyr::replace_na('BRCA1')
combined_results$variant_type  = combined_results$variant_type %>% tidyr::replace_na('frameshift')
combined_results$HGVS_united  = combined_results$HGVS_united %>% tidyr::replace_na('p.Glu23ValfsTer17')

print(combined_results, width=Inf)
#quit()

# compare with clinical record ----
britroc_con = make_connection_to_postgres_server('britroc1','jblab-db.cri.camres.org',5432)
patient_table = dbReadTable(britroc_con, 'patients') %>% tibble::as_tibble() 

clinical_hrd_patients = patient_table %>% 
	dplyr::filter(brca_status %in% c('BRCA1 mutation identified','BRCA2 mutation identified')) %>%
	.$britroc_number	
tamseq_hrd_patients = combined_results %>% dplyr::filter(final_call_set==TRUE) %>% .$fk_britroc_number %>% unique()

print(clinical_hrd_patients)
print(tamseq_hrd_patients)

clinical_hrd_patients_exclusive = setdiff(clinical_hrd_patients,tamseq_hrd_patients)
patients_with_no_germline_sequencing = 
	setdiff(1:276, readr::read_tsv('config/germline_metadata.tsv') %>% .$fk_britroc_number)
#print(patients_with_no_germline_sequencing)

clinical_hrd_patients_exclusive = setdiff(clinical_hrd_patients_exclusive, patients_with_no_germline_sequencing)

print(clinical_hrd_patients_exclusive)
 
intersect(clinical_hrd_patients, patients_with_no_germline_sequencing)

# add clinical records to the existing table ----
colnames(combined_results) %>% print()

new_combined_results_records = tibble::tribble(
	 ~fk_britroc_number, ~gene_symbol, ~chromosome, ~position, ~ref, ~alt, ~variant_type, ~HGVS_united, ~called_by_HaplotypeCaller, ~called_by_octopus, ~germline_sequencing, ~final_call_set, ~comment,
	46, 'BRCA1', 'chr17', NA, NA, NA, NA, NA, FALSE, FALSE, TRUE, TRUE, 'No additional variant information provided from clinical records',
	86, 'BRCA2', 'chr13', 32912654, 'CT', 'A', 'frameshift', 'p.Thr1388Asnfs*22', FALSE, FALSE, TRUE, TRUE, 'Evidence from clinical records; some inconclusive evidence from IGV inspection',
	127, 'BRCA1', 'chr17', NA, NA, NA, 'frameshift', 'p.Ala145_Asp1662fsX14', FALSE, FALSE, TRUE, TRUE, 'Multi-exon deletion reported in clinical records; impossible for thisto be detected by amplicon sequencing',
	175, 'BRCA2', 'chr13', NA, NA, NA, NA, NA, FALSE, FALSE, TRUE, TRUE, 'Multi-exon deletion reported in clinical records; impossible for this to be detected by amplicon sequencing',
	225, 'BRCA1', 'chr17', NA, NA, NA, NA, NA, NA, NA, FALSE, TRUE,  'Evidence from clinical records only',
	245, 'BRCA2', 'chr13', NA, NA, NA, NA, NA, FALSE, FALSE, TRUE, TRUE, 'No additional variant information provided from clinical records',
	247, 'BRCA1', 'chr17', NA, NA, NA, NA, NA, NA, NA, FALSE, TRUE, 'Evidence from clinical records only',
	257, 'BRCA2', 'chr13', NA, NA, NA, NA, NA, NA, NA, FALSE, TRUE, 'Evidence from clinical records only',
	259, 'BRCA1', 'chr17', NA, NA, NA, NA, NA, NA, NA, FALSE, TRUE, 'Evidence from clinical records only',
	275, 'BRCA1', 'chr17', NA, NA, NA, NA, NA, NA, NA, FALSE, TRUE, 'Evidence from clinical records only',
	276, 'BRCA2', 'chr13', NA, NA, NA, NA, NA, NA, NA, FALSE, TRUE, 'Evidence from clinical records only' 
)

print(new_combined_results_records, width=Inf)

combined_results = rbind(combined_results, new_combined_results_records)

combined_results = combined_results %>% dplyr::arrange(-final_call_set, fk_britroc_number)

print(combined_results)

readr::write_tsv(combined_results, 'BriTROC-1_germline_variants.tsv', col_names=TRUE)

combined_results$fk_britroc_number %>% table %>% sort()

combined_results %>% dplyr::filter(final_call_set==TRUE) %>% .$HGVS_united %>% table() %>% sort()

quit()





