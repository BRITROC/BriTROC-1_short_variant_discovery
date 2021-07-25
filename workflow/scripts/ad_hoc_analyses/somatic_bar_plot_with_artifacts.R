
library(magrittr)
library(ggplot2)
library(patchwork)

somatic_mutations = tibble::tribble(
	~Gene, ~num_patients, ~tumour_type,
	'RAD51C',     0L, 'archival',
	'RAD51C',	0L, 'relapse',
	'RAD51D',     0L, 'archival',
	'RAD51D',     1L, 'relapse',
	'BARD1',      2L, 'archival',
	'BARD1',	1L, 'relapse', 
	'RAD51B',     0L, 'archival',
	'RAD51B', 	0L, 'relapse',
	'BRIP1',      0L, 'archival',
	'BRIP1',	2L, 'relapse',
	'PALB2',      0L, 'archival',
	'PALB2',	2L, 'relapse', 
	'BRCA1',     1L, 'archival',
	'BRCA1',	1L, 'relapse', 
	'FANCM',     3L, 'archival',
	'FANCM',	2L, 'relapse',
	'BRCA2',     4L, 'archival',
	'BRCA2',	3L, 'relapse'
)

somatic_mutations$prop_patients = somatic_mutations$num_patients / 141

num_patients = 141

p1 = ggplot(
	somatic_mutations, aes(x=Gene, y=num_patients,fill=tumour_type)
) + geom_col(position = position_dodge2(preserve = "single")) + 
   ylab('Number of patients') + 
   ggtitle('BriTROC-1: Potentially deleterious short somatic mutations') +
   labs(subtitle='141 patients')

p2 = ggplot(
        somatic_mutations, aes(x=Gene, y=prop_patients,fill=tumour_type)
) + geom_col(position = position_dodge2(preserve = "single")) + 
	ylab('Proportion of patients') + 
	ggtitle('BriTROC-1: Potentially deleterious short somatic mutations') + 
	labs(subtitle='141 patients')

p3 = p1 + p2

ggsave(p3, filename='britroc1_somatic_mutations.png', device='png', scale=1.0, width=16)

print(somatic_mutations)
