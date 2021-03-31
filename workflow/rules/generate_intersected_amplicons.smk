rule all:
	input: 'resources/intersected_panel_6_28_amplicons.txt'

rule get_intersected_amplicon_file:
	input:
		panel_6='/scratcha/jblab/amplicon_panels/6_HR_Elkes_panel_2015_08_28/targets.txt',
		panel_28='/scratcha/jblab/amplicon_panels/28_JBLAB_AAprimers_dream_panel/targets.txt'
	output: 'resources/intersected_panel_6_28_amplicons.interval_list'
	script: '../scripts/intersect_amplicons.R' 
