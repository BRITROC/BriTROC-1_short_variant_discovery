# generate amplicon and targets intervals file which is the intersection of amplicons and targets for panel 6 and panel 28

rule intersect_intervals:
	input:
		panel_6='resources/amplicon_panels/6_HR_Elkes_panel_2015_08_28/{interval_type}.txt',
		panel_28='resources/amplicon_panels/28_JBLAB_AAprimers_dream_panel/{interval_type}.txt'
	output: temp('resources/intersected_panel_6_28_amplicons.{interval_type}.interval_list.tmp')
	script: '../scripts/generate_interected_amplicons/intersect_amplicons.R'

rule antijoin_intervals:
	input:
		panel_6='resources/amplicon_panels/6_HR_Elkes_panel_2015_08_28/{interval_type}.txt',
		panel_28='resources/amplicon_panels/28_JBLAB_AAprimers_dream_panel/{interval_type}.txt'
	output: temp('resources/antijoined_panel_28_6_amplicons.{interval_type}.interval_list.tmp')
	script: '../scripts/generate_intersected_amplicons/antijoin_amplicons.R'

rule union_of_tp53_intervals:
	input:
		panel_1='resources/amplicon_panels/1_JBLAB_AAprimers_48_basic_panel/{interval_type}.txt',
		panel_10='resources/amplicon_panels/10_JBLAB_AAprimers_TP53short_only/{interval_type}.txt',
		panel_28='resources/amplicon_panels/28_JBLAB_AAprimers_dream_panel/{interval_type}.txt',
	output: temp('resources/union_of_tp53_amplicons.{interval_type}.interval_list.tmp')
	script: '../scripts/generate_intersected_amplicons/get_union_of_tp53_intervals.R'

rule add_header_to_interval_list_file:
	input:
		template='resources/amplicon_panels/28_JBLAB_AAprimers_dream_panel/{interval_type}.txt',
		incomplete_intervals_file='resources/{amplicon_panel}.{interval_type}.interval_list.tmp'
	output: 'resources/{amplicon_panel}.{interval_type}.interval_list'
	shell: "grep '@' {input.template} | cat - {input.incomplete_intervals_file} > temp_file && mv temp_file {output}"
