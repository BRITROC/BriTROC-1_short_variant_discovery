# generate amplicon and targets intervals file which is the intersection of amplicons and targets for panel 6 and panel 28

rule intersect_intervals:
	input:
		panel_6='/scratcha/jblab/amplicon_panels/6_HR_Elkes_panel_2015_08_28/{interval_type}.txt',
		panel_28='/scratcha/jblab/amplicon_panels/28_JBLAB_AAprimers_dream_panel/{interval_type}.txt'
	output: temp('resources/intersected_panel_6_28_amplicons.{interval_type}.interval_list.tmp')
	script: '../scripts/intersect_amplicons.R'

rule add_header_to_interval_list_file:
	input: 
		template='/scratcha/jblab/amplicon_panels/28_JBLAB_AAprimers_dream_panel/{interval_type}.txt',
		incomplete_intervals_file=rules.intersect_intervals.output
	output: 'resources/intersected_panel_6_28_amplicons.{interval_type}.interval_list'
	shell: "grep '@' {input.template} | cat - {input.incomplete_intervals_file} > temp_file && mv temp_file {output}"
