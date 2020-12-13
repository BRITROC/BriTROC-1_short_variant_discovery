rule all:
	input: 'intersected_panel_6_28_amplicons.txt'

rule get_intersected_amplicon_file:
	input:
		panel_6='/scratcha/jblab/amplicon_panels/6_HR_Elkes_panel_2015_08_28/targets.txt',
		panel_28='/scratcha/jblab/amplicon_panels/28_JBLAB_AAprimers_dream_panel/targets.txt'
	output: 'intersected_panel_6_28_amplicons.interval_list'
	script: 'intersect_amplicons.R' 
		
		#"awk 'NR==FNR{{A[$1];next}}$5 in A' {input.panel_6} {input.panel_28} > {output}"
